import os
import numpy as np

import frb_olympics


class bonsai_dedisperser(frb_olympics.dedisperser_base):
    def __init__(self, config, name=None):
        """
        The 'config' argument can be either a dictionary or a filename.
        """

        try:
            import bonsai
        except:
            raise ImportError("frb_olympics: couldn't import 'bonsai'.  You may need to install it from https://github.com/CHIMEFRB/ch_frb_io.")

        if name is None:
            name = config if isinstance(config, basename) else 'bonsai_dedisperser'
            name = os.path.basename(name)
            name = name.split('.')[0]
        
        frb_olympics.dedisperser_base.__init__(name)

        self.dedisperser = bonsai.Dedisperser(config, fill_rfi_mask=False, allocate=False, use_analytic_normalization=True)


    def init_search_params(self, search_params):
        if hasattr(self, 'search_params'):
            raise RuntimeError('double call to frb_olympics.bonsai_dedisperser.init_search_params()')

        # Lots of sanity checks
        assert isinstance(sp, frb_olympics.search_params)
        assert sp.nfreq == self.dedisperser.nfreq
        assert sp.dm_max <= np.max(self.dedisperser.max_dm)
        assert abs(sp.freq_lo_MHz - self.dedisperser.freq_lo_MHz) < 1.0e-3
        assert abs(sp.freq_hi_MHz - self.dedisperser.freq_hi_MHz) < 1.0e-3
        assert abs(sp.dt_sample - self.dedisperser.dt_sample) < 1.0e-7

        import bonsai

        self.search_params = sp
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
                 'tfin': self.global_max_tracker.global_max_trigger_arrival_time }


    def deallocate(self):
        self.dedisperser.deallocate()


    def jsonize(self):
        return self.dedisperser.config


    @staticmethod
    def from_json(j, filename=None):
        (j, filename) = frb_olympics._from_json_helper(j, filename, 'bonsai_dedisperser.from_json()')

        j = copy.copy(j)
        j.pop('module_name')
        j.pop('class_name')

        name = j.pop('name', filename)
        return bonsai_dedispeser(j, name)
