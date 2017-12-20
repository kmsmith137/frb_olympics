import os
import copy
import numpy as np

import frb_olympics


class bonsai_dedisperser(frb_olympics.dedisperser_base):
    def __init__(self, config, tex_label):
        """
        The 'config' argument can be either a dictionary or a filename.
        """

        try:
            import bonsai
        except:
            raise ImportError("frb_olympics: couldn't import 'bonsai'.  You may need to install it from https://github.com/CHIMEFRB/ch_frb_io.")

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

            intensity[:,:(it1-it0)] = arr[:,it0:it1]
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

        config = copy.copy(r.json)
        config.pop('module_name')
        config.pop('class_name')

        tex_label = config.pop('tex_label')
        return bonsai_dedisperser(config, tex_label)
