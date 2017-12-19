import os
import numpy as np

from .frb_olympics import search_params, dedisperser_base


class bonsai_dedisperser(dedisperser_base):
    def __init__(self, config_filename):
        try:
            import bonsai
        except:
            raise ImportError("frb_olympics: couldn't import 'bonsai'.  You may need to install it from https://github.com/CHIMEFRB/ch_frb_io.")

        self.config_filename = os.path.abspath(config_filename)
        
        self.dedisperser = bonsai.Dedisperser(config_filename,
                                              fill_rfi_mask = False,
                                              allocate = False,
                                              use_analytic_normalization = True)


    def init_search_params(self, sp):
        assert isinstance(sp, search_params)
        self.search_params = sp

        # Lots of sanity checks
        assert sp.nfreq == self.dedisperser.nfreq
        assert sp.dm_max <= np.max(self.dedisperser.max_dm)
        assert abs(sp.freq_lo_MHz - self.dedisperser.freq_lo_MHz) < 1.0e-3
        assert abs(sp.freq_hi_MHz - self.dedisperser.freq_hi_MHz) < 1.0e-3
        assert abs(sp.dt_sample - self.dedisperser.dt_sample) < 1.0e-7

        import bonsai
        self.global_max_tracker = bonsai.global_max_tracker(sp.dm_min, sp.dm_max)
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
                 'tf': self.global_max_tracker.global_max_trigger_arrival_time }


    def deallocate(self):
        self.dedisperser.deallocate()


    def jsonize(self):
        return {
            'module_name': self.__module__,
            'class_name': self.__class__.__name__,
            'config_filename': self.config_filename
        }


    @staticmethod
    def from_json(j, json_filename):
        return bonsai_dedisperser(j['config_filename'])
