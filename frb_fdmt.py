import cpuFDMT
import numpy as np


def is_power_of_two(n):
    return (n >= 1) and ((n & (n-1)) == 0)


class fdmt:
    def __init__(self, maxDT, recdepth, doWeighting=True):
        # The search algorithm onstructor should initialize self.name, in addition 
        # to any fields needed internally.

        # I wasn't sure if FDMT had the same power-of-two restriction that the tree
        # search has, so I added this assert to be on the safe side...
        assert is_power_of_two(maxDT)
        assert recdepth >= 0

        self.name = ('fdmt-%d-rec%d' % (maxDT,recdepth))
        self.maxDT = maxDT
        self.recdepth = recdepth
        self.doWeighting = doWeighting

        if not doWeighting:
            self.name += '-noweight'


    def search_init(self, search_params):
        # search_init() is called when the search_params are specified.
        # It should initialize: self.search_gb, self.debug_buffer_ndm, self.debug_buffer_nt
        # (See docstring at top of frb_olympics.py)

        # Incremental FDMT search not supported yet, so we just raise an exception
        # if search_params.nchunks is larger than 1 (indicating an incremental search)
        assert search_params.nchunks == 1

        self.search_params = search_params
        self.search_gb = 1.0e-9 * self.maxDT * search_params.nsamples_per_chunk * 4.0   # approximate
        self.search_result = 0.0

        self.debug_buffer_ndm = self.maxDT
        self.debug_buffer_nt = search_params.nsamples_per_chunk // 2**self.recdepth

        search_max_delay = search_params.dt_sample * (self.maxDT-1) * 2**self.recdepth   # seconds
        search_max_dm = search_max_delay / 4.148806e3 / (search_params.band_freq_lo_MHz**(-2) - search_params.band_freq_hi_MHz**(-2))
        print '%s: initialized, max DM requested = %s, max DM searched by FDMT = %s' % (self.name, search_params.dm_max, search_max_dm)

        # If the FDMT hasn't been configured to search the full DM range of the sims, 
        # then the result of the frb_olympics won't be very meaningful, so just abort.
        assert search_params.dm_max <= search_max_dm 


    def search_start(self, mpi_rank_within_node):
        #
        # search_start() is called just before the timestream is processed
        # with search_chunk().  It usually just allocates some buffers which will be
        # deallocated later in search_end().  (See docstring at top of frb_olympics.py)
        #
        # In the case of the FDMT, we use search_start() to initialize some fields
        # which are globals, rather than arguments to cpuFDMT.fdmt().
        #
        cpuFDMT.fmin = self.search_params.band_freq_lo_MHz
        cpuFDMT.fmax = self.search_params.band_freq_hi_MHz
        cpuFDMT.nchan = self.search_params.nchan
        cpuFDMT.maxDT = self.maxDT
        cpuFDMT.doWeighting = self.doWeighting
        cpuFDMT.fs, cpuFDMT.df = np.linspace(cpuFDMT.fmin, cpuFDMT.fmax, cpuFDMT.nchan, endpoint=False, retstep=True)
        cpuFDMT.A = None
        cpuFDMT.B = None
        cpuFDMT.Q = None


    def search_chunk(self, chunk, ichunk, debug_buffer=None):
        #
        # search_chunk() is called to search a timestream chunk for an FRB.  It sets 
        # self.search_result to the "trigger" statistic T from the FRB olympics memo.  
        #
        # In a non-incremental search, search_chunk() will be called once with ichunk=0.
        #
        # In an incremental search, search_chunk() will be called multiple times (with 
        # ichunk=0,..,nchunks-1) and search_chunk() is responsible for saving state between
        # calls and setting self.search_result at the end.
        #
        # Normally, search_chunk() just computes the trigger statistic (a scalar) without
        # actually returning the output of the DM transform (a 2D array indexed by DM and
        # arrival time).  However, if debug_buffer is not None, it will be an array of shape
        # (self.debug_buffer_ndm, self.debug_buffer_nt) which holds the output of the DM
        # transform.  I usually use this for visual inspection with frb-dump.py, but it
        # could be used for other things.
        #
        # "Masked" entries in the debug_buffer which are not actually searched
        # (because their DM or arrival time is out of range) are indicated by setting
        # them to -1.0e30.
        #

        # incremental search not supported yet
        # (this assert is redundant, since the same assert appears in search_init(), but that's ok!)
        assert self.search_params.nchunks == 1 and ichunk == 0

        assert chunk.shape == (self.search_params.nchan, self.search_params.nsamples_per_chunk)

        if debug_buffer is not None:
            #
            # The FDMT search runs a few times at different levels of downsampling (see 
            # cpuFDMT.recursive_fdmt()).  It made the most sense to me to fill the debug_buffer
            # with the output of the DM transform at the highest level of downsampling.  This
            # is because the searched region in the (DM,arrival_time) space is largest, so it
            # makes the most informative plot for visual debugging purposes.
            #
            I = chunk
            for i in xrange(self.recdepth):
                # cut-and-paste from cpuFDMT.py
                I = I[:,::2]+I[:,1::2] if (I.shape[1]%2 == 0) else I[:,:-1:2]+I[:,1::2]
                if self.doWeighting: 
                    I /= np.sqrt(2.)
            
            assert I.shape == (self.search_params.nchan, self.debug_buffer_nt)

            # FDMT ends up being run twice, which is inefficient, but that's OK on a debug path
            debug_buffer[:,:] = cpuFDMT.fdmt(I, returnMaxSigma=False)[:,:]
            
            # mask entries which aren't searched by fdmt(), for consistency with cpuFDMT.fdmt()
            for i in xrange(debug_buffer.shape[0]):
                debug_buffer[i,:i] = -1.0e30

        self.search_result = cpuFDMT.recursive_fdmt(chunk, depth=self.recdepth)


    def search_end(self):
        #
        # search_end() just deallocates buffers.
        #
        # This is important so that we don't use too much memory when multiple searches are
        # run in parallel with frb-compare.py.
        #
        cpuFDMT.A = None
        cpuFDMT.B = None
        cpuFDMT.Q = None
