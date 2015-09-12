import cpuFDMT
import numpy as np


def is_power_of_two(n):
    return (n >= 1) and ((n & (n-1)) == 0)


class fdmt:
    def __init__(self, maxDT, recdepth, doWeighting=True):
        assert is_power_of_two(maxDT)  # not sure if this is necessary
        assert recdepth >= 0

        self.name = ('fdmt-%d-rec%d' % (maxDT,recdepth))
        self.maxDT = maxDT
        self.recdepth = recdepth
        self.doWeighting = doWeighting

        if not doWeighting:
            self.name += '-noweight'


    def search_init(self, search_params):
        # incremental search not implemented yet
        assert search_params.nchunks == 1

        self.search_params = search_params
        self.search_gb = 1.0e-9 * self.maxDT * search_params.nsamples_per_chunk * 4.0   # approximate
        self.search_result = 0.0

        self.debug_buffer_ndm = self.maxDT
        self.debug_buffer_nt = search_params.nsamples_per_chunk // 2**self.recdepth

        search_max_delay = search_params.dt_sample * (self.maxDT-1) * 2**self.recdepth   # seconds
        search_max_dm = search_max_delay / 4.148806e3 / (search_params.band_freq_lo_MHz**(-2) - search_params.band_freq_hi_MHz**(-2))
        print '%s: initialized, max DM requested = %s, max DM searched by FDMT = %s' % (self.name, search_params.dm_max, search_max_dm)

        # if the FDMT hasn't been configured to search the full DM range of the sims, 
        # then the result of the frb_olympics won't be very meaningful, so just abort
        assert search_params.dm_max <= search_max_dm 


    def search_start(self):
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
        # incremental search not supported yet
        assert self.search_params.nchunks == 1

        if debug_buffer is not None:
            # FDMT ends up being run twice, which is inefficient, but that's OK on a debug path
            I = chunk
            for i in xrange(self.recdepth):
                # cut-and-paste from cpuFDMT.py
                I = I[:,::2]+I[:,1::2] if (I.shape[1]%2 == 0) else I[:,:-1:2]+I[:,1::2]
                if self.doWeighting: 
                    I /= np.sqrt(2.)
            
            assert I.shape == (self.search_params.nchan, self.debug_buffer_nt)
            debug_buffer[:,:] = cpuFDMT.fdmt(I, returnMaxSigma=False)[:,:]
            
            # mask entries which aren't searched by fdmt()
            for i in xrange(debug_buffer.shape[0]):
                debug_buffer[i,:i] = -1.0e30

        self.search_result = cpuFDMT.recursive_fdmt(chunk, depth=self.recdepth)


    def search_end(self):
        cpuFDMT.A = None
        cpuFDMT.B = None
        cpuFDMT.Q = None
