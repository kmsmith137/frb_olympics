import numpy as np
from frb_olympics_c import frb_search_params


class downsample:
    def __init__(self, algo, ndownsample):
        assert ndownsample >= 1
        self.algo = algo
        self.nds = ndownsample
        self.name = '%s-ds%d' % (algo.name, ndownsample)


    def search_init(self, search_params):
        assert search_params.nsamples_per_chunk % self.nds == 0

        # copy and modify
        p = frb_search_params(search_params)
        p.dt_sample = p.dt_sample * self.nds
        p.nsamples_tot = p.nsamples_tot // self.nds
        p.nsamples_per_chunk = p.nsamples_per_chunk // self.nds

        self.algo.search_init(p)
        self.search_params = search_params
        self.debug_buffer_ndm = self.algo.debug_buffer_ndm
        self.debug_buffer_nt = self.algo.debug_buffer_nt
        self.search_gb = self.algo.search_gb + (1.0e-9 * p.nchan * p.nsamples_per_chunk * 4)


    def search_start(self):
        self.algo.search_start()
        self.search_result = self.algo.search_result


    def search_chunk(self, chunk, ichunk, debug_buffer=None):
        nchan = self.search_params.nchan
        nt = self.search_params.nsamples_per_chunk // self.nds
        assert chunk.shape == (nchan, self.nds * nt)

        t = np.zeros((nchan,nt), dtype=np.float32)
        for i in xrange(self.nds):
            t[:,:] += chunk[:,i::self.nds]

        t[:,:] /= float(self.nds)**0.5
        self.algo.search_chunk(t, ichunk, debug_buffer)
        self.search_result = self.algo.search_result


    def search_end(self):
        self.algo.search_end()
        self.search_result = self.algo.search_result
