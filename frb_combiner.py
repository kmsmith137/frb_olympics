import numpy as np
from frb_olympics_c import frb_search_params


class combiner:
    """
    A wrapper algorithm which combines multiple "child" algorithms,
    and takes the max of the search_results.
    """

    def __init__(self, name=None):
        self.name = name if name is not None else 'combine'

        # list of (algo, search_params_kwds) pairs
        self.algo_list = [ ]


    def add_algo(self, algo, **kwds):
        self.algo_list.append((algo, kwds))


    def search_init(self, search_params):
        assert len(self.algo_list) >= 1

        for (algo, kwds) in self.algo_list:
            # copy and modify
            p = frb_search_params(search_params)
            for (k,v) in kwds.iteritems():
                setattr(p, k, v)

            algo.search_init(p)

        algo0 = self.algo_list[0][0]
        self.search_params = search_params
        self.debug_buffer_ndm = algo0.debug_buffer_ndm
        self.debug_buffer_nt = algo0.debug_buffer_nt
        self.search_gb = sum(algo.search_gb for (algo,kwds) in self.algo_list)


    def search_start(self, mpi_rank_within_node):
        for (algo, kwds) in self.algo_list:
            algo.search_start(mpi_rank_within_node)

        self.search_result = max(algo.search_result for (algo,kwds) in self.algo_list)


    def search_chunk(self, chunk, ichunk, debug_buffer=None):
        for (i,(algo,kwds)) in enumerate(self.algo_list):
            d = debug_buffer if (i == 0) else None
            algo.search_chunk(chunk, ichunk, debug_buffer)

        self.search_result = max(algo.search_result for (algo,kwds) in self.algo_list)        


    def search_end(self):
        for (algo,kwds) in self.algo_list:
            algo.search_end()

        self.search_result = max(algo.search_result for (algo,kwds) in self.algo_list)

