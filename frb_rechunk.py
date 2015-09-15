import numpy as np
from frb_olympics_c import frb_search_params


class rechunk:
    """
    A wrapper algorithm which rechunks the data using a different chunk size.

    I usually use this to debug incremental search, by verifying that the
    same algorithm run with different chunk sizes gives the same answer.
    """

    def __init__(self, algo, new_nsamples_per_chunk):
        assert new_nsamples_per_chunk > 0

        self.algo = algo
        self.new_nsamples_per_chunk = new_nsamples_per_chunk
        self.name = '%s-rechunk%d' % (algo.name, new_nsamples_per_chunk)


    def search_init(self, search_params):
        assert search_params.nsamples_tot % self.new_nsamples_per_chunk == 0

        # copy and modify
        p = frb_search_params(search_params)
        p.nsamples_per_chunk = self.new_nsamples_per_chunk
        p.nchunks = p.nsamples_tot // self.new_nsamples_per_chunk

        self.algo.search_init(p)
        self.search_params = search_params
        self.debug_buffer_ndm = self.algo.debug_buffer_ndm
        self.debug_buffer_nt = self.algo.debug_buffer_nt
        self.search_gb = self.algo.search_gb + (1.0e-9 * p.nchan * self.new_nsamples_per_chunk * 4)


    def search_start(self):
        self.algo.search_start()
        self.search_result = self.algo.search_result
        
        self.tbuf = np.zeros((self.search_params.nchan, self.new_nsamples_per_chunk), dtype=np.float32)
        self.tbuf_curr = 0   # number of timeslices currently in buffer
        self.tcount = 0      # number of times self.algo.search_chunk() has been called


    def search_chunk(self, chunk, ichunk, debug_buffer=None):
        nt_in = self.search_params.nsamples_per_chunk
        nt_out = self.new_nsamples_per_chunk
        assert chunk.shape == (self.search_params.nchan, nt_in)

        s = 0   # current position in source buffer
        while s < nt_in:
            d = self.tbuf_curr    # current position in dest buffer
            n = min(nt_in - s, nt_out - d)
            self.tbuf[:,d:d+n] = chunk[:,s:s+n]
            self.tbuf_curr += n
            s += n

            if self.tbuf_curr == nt_out:
                self.algo.search_chunk(self.tbuf, self.tcount, debug_buffer)
                self.search_result = self.algo.search_result
                self.tbuf_curr = 0
                self.tcount += 1


    def search_end(self):
        assert self.tcount == self.search_params.nsamples_tot // self.new_nsamples_per_chunk
        self.algo.search_end()
        self.search_result = self.algo.search_result
