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
        The 'config' argument can be either a dictionary or a filename.
        """

        if not import_successful:
            # Rather than throw an exception, we let 'import bonsai' throw an uncaught
            # exception, so that the caller can see what the problem is.
            import bonsai
            raise RuntimeError("frb_olympics.bonsai_dedisperser internal error: 'import bonsai' worked on the second try?!")

        frb_olympics.dedisperser_base.__init__(self, tex_label)

        self.dedisperser = bonsai.Dedisperser(config, fill_rfi_mask=False, allocate=False, use_analytic_normalization=True)


    def init_search_params(self, sparams):
        if hasattr(self, 'search_params'):
            raise RuntimeError('double call to frb_olympics.bonsai_dedisperser.init_search_params()')

        # Lots of sanity checks
        assert isinstance(sparams, frb_olympics.search_params)
        assert sparams.nfreq == self.dedisperser.nfreq
        assert sparams.dm_max <= np.max(self.dedisperser.max_dm)
        assert abs(sparams.freq_lo_MHz - self.dedisperser.freq_lo_MHz) < 1.0e-3
        assert abs(sparams.freq_hi_MHz - self.dedisperser.freq_hi_MHz) < 1.0e-3
        assert abs(sparams.dt_sample - self.dedisperser.dt_sample) < 1.0e-7

        import bonsai

        self.search_params = sparams
        self.global_max_tracker = bonsai.global_max_tracker(sparams.dm_min, sparams.dm_max)
        self.dedisperser.add_processor(self.global_max_tracker)


    def allocate(self):
        self.dedisperser.allocate()


    def dedisperse(self, arr):
        nfreq = self.search_params.nfreq
        nsamples = self.search_params.nsamples
        nt_chunk = self.dedisperser.nt_chunk

        assert arr.shape == (nfreq, nsamples)

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
        self.dedisperser.deallocate()


    def jsonize(self):
        return self.dedisperser.config


    @staticmethod
    def from_json(j, filename=None):
        r = frb_olympics.json_read_helper(j, filename, 'bonsai_dedipserser.from_json()')
        j = copy.copy(r.json)

        j.pop('module_name')
        j.pop('class_name')
        tex_label = j.pop('tex_label')

        return bonsai_dedisperser(j, tex_label)
