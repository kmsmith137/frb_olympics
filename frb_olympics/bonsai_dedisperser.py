"""
bonsai dedisperser (https://github.com/CHIMEFRB/bonsai).

Note that the DM's and arrival times reported by bonsai are usually coarse-grained for speed, and 
therefore contain errors due to coarse-graining.  To disable the coarse-graining, you can set
`dm_coarse_graining` and `time_coarse_graining` to 1 in the bonsai config file.

To do list:

  - Currently, the analytic transfer matrix is computed from scratch whenever a bonsai_dedisperser
    is constructed, which is annoying!  Should cache it in an HDF5 file (this is mostly implemented
    in bonsai already.)

  - Cleanup: less verbose output from bonsai_dedisperser.jsonize()
"""


import os
import copy
import numpy as np

import frb_olympics


####################################################################################################


import_successful = False

try:
    import bonsai
    import_successful = True
except ImportError:
    pass


####################################################################################################


class bonsai_dedisperser(frb_olympics.dedisperser_base):
    def __init__(self, config, tex_label):
        """
        The 'config' argument can be either a dictionary containing the bonsai config_params
        (there are a lot of these, see bonsai documentation for more info!), or the filename
        of a valid bonsai config file.
        """

        if not import_successful:
            # Rather than throw an exception, we let 'import bonsai' throw an uncaught
            # exception, so that the caller can see what the problem is.

            import bonsai as b
            raise RuntimeError("frb_olympics.bonsai_dedisperser internal error: 'import bonsai' worked on the second try?!")

        # Calling the bonsai.Dedisperser constructor with use_analytic_normalization=True 
        # means that bonsai will precompute the exact variance of its output array, assuming 
        # that each element of the input intensity array is a unit Gaussian.  The bonsai
        # output array is then normalized to "sigmas".
        #
        # This setting is unsuitable for real data, but well-suited for the frb_olympics,
        # where we know that the noise is just an array of unit Gaussians.  In particular
        # it means that it is safe to use 'run-frb-olympics' with the -N flag, which removes
        # noise from the simulations (for more info, see 'run-frb-olympics -h').
        #
        # Note that because we're using the analytic normalization, we can call the
        # dedipserser_base subclass constructor with precomputed_variance=True.

        frb_olympics.dedisperser_base.__init__(self, tex_label, precomputed_variance=True)

        self.dedisperser = bonsai.Dedisperser(config, fill_rfi_mask=False, allocate=False, use_analytic_normalization=True)


    def init_search_params(self, sparams):
        """Overrides dedisperser_base.init_search_params().  The 'sparams' argument is an instance of 'class search_params'."""

        if hasattr(self, 'search_params'):
            raise RuntimeError('double call to frb_olympics.bonsai_dedisperser.init_search_params()')

        # Lots of sanity checks
        assert isinstance(sparams, frb_olympics.search_params)
        assert sparams.nfreq == self.dedisperser.nfreq
        assert sparams.dm_max <= np.max(self.dedisperser.max_dm)
        assert abs(sparams.freq_lo_MHz - self.dedisperser.freq_lo_MHz) < 1.0e-3
        assert abs(sparams.freq_hi_MHz - self.dedisperser.freq_hi_MHz) < 1.0e-3
        assert abs(sparams.dt_sample - self.dedisperser.dt_sample) < 1.0e-7

        self.search_params = sparams
        self.global_max_tracker = bonsai.global_max_tracker(sparams.dm_min, sparams.dm_max)
        self.dedisperser.add_processor(self.global_max_tracker)


    def allocate(self):
        """Overrides dedisperser_base.allocate()."""
        self.dedisperser.allocate()


    def dedisperse(self, arr):
        """
        Overrides dedisperser_base.dedisperse().  The 'arr' argument is a float32 array of shape (nfreq, nsamples).

        The return value of dedisperse() is a dictionary with 3 members: 'snr', 'dm', and one of { 'tmid', 'tini', 'tfin' }.
        This allows the dedisperser to use one of three possible definitions of the arrival time of an FRB:
           - "initial" arrival time: arrival time at the highest frequency in the band (i.e. least delayed)
           - "final" arrival time: arrival time at the lowest frequency in the band (i.e. most delayed)
           - "middle" arrival time: average of initial and final times (warning: not the arrival time at the central frequency!)

        We return 'tfin', since this is what bonsai "naturally" computes.
        """

        nfreq = self.search_params.nfreq
        nsamples = self.search_params.nsamples
        nt_chunk = self.dedisperser.nt_chunk

        assert arr.shape == (nfreq, nsamples)

        # Note that bonsai is designed to operate incrementally in chunks of fixed length 'nt_chunk', 
        # so we divide the timestream into chunks.  If the timestream ends in the middle of a chunk,
        # we indicate this by setting the bonsai weights array to zero.

        for it0 in xrange(0, nsamples, nt_chunk):
            it1 = min(it0 + nt_chunk, nsamples)

            intensity = np.zeros((nfreq, nt_chunk), dtype=np.float32)
            weights = np.zeros((nfreq, nt_chunk), dtype=np.float32)

            # Note that the frequency channel ordering gets reversed here: frb_olympics uses
            # lowest-to-highest channel ordering, and bonsai uses the opposite.

            intensity[:,:(it1-it0)] = arr[::-1,it0:it1]
            weights[:,:(it1-it0)] = 1.0
            
            self.dedisperser.run(intensity, weights)

        self.dedisperser.end_dedispersion()

        return { 'snr': self.global_max_tracker.global_max_trigger,
                 'dm': self.global_max_tracker.global_max_trigger_dm,
                 'tfin': self.global_max_tracker.global_max_trigger_arrival_time }


    def deallocate(self):
        """Overrides dedisperser_base.deallocate()."""
        self.dedisperser.deallocate()


    def jsonize(self):
        """Overrides dedisperser_base.jsonize().  The return value should be a python dictionary which is valid JSON."""
        return self.dedisperser.config


    @staticmethod
    def from_json(j, filename=None):
        """
        This is the "inverse" of jsonize(): it takes a python dictionary which is valid JSON, 
        and returns a bonsai_dedisperser instance.
        """

        r = frb_olympics.json_read_helper(j, filename, 'bonsai_dedipserser.from_json()')
        j = copy.copy(r.json)

        # Detail: These keys are added by the caller (i.e. they are not part of the return value of
        # bonsai_dedisperser.jsonize()).  We remove them from the dictionary to avoid getting warnings
        # from bonsai.

        j.pop('module_name')
        j.pop('class_name')
        tex_label = j.pop('tex_label')

        return bonsai_dedisperser(j, tex_label)
