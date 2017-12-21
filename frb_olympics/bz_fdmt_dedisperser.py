import sys
import numpy as np
import frb_olympics


####################################################################################################


import_successful = False

try:
    import bz_fdmt
    import_successful = True
except ImportError:
    pass


####################################################################################################


def is_power_of_two(n):
    return (n > 0) and (n & (n-1)) == 0

def round_up_to_power_of_two(n):
    assert n > 0
    return 2**(int(np.log2(n-0.5)) + 1)


class bz_fdmt_dedisperser(frb_olympics.dedisperser_base):
    def __init__(self):
        if not import_successful:
            # Rather than throw an exception, we let 'import bz_fdmt' throw an uncaught
            # exception, so that the caller can see what the problem is.
            import bz_fdmt
            raise RuntimeError("frb_olympics.bz_fdmt_dedisperser internal error: 'import bz_fdmt' worked on the second try?!")

        frb_olympics.dedisperser_base.__init__(self, tex_label='FDMT')


    def init_search_params(self, sparams):
        self.search_params = sparams
        self.freq_lo = sparams.freq_lo_MHz
        self.freq_hi = sparams.freq_hi_MHz
        self.dt_sample = sparams.dt_sample

        # These restrictions can be relaxed.
        if not is_power_of_two(sparams.nfreq):
            raise RuntimeError("frb_olympics.bz_fdmt_dedisperser: we currently require nfreq to be a power of two")
        if not is_power_of_two(sparams.nsamples):
            raise RuntimeError("frb_olympics.bz_fdmt_dedisperser: we currently require nt_in to be a power of two")

        # samples_per_dm = dispersion delay at DM=1, in samples
        dt1_i = frb_olympics.dispersion_delay(1.0, self.freq_hi)
        dt1_f = frb_olympics.dispersion_delay(1.0, self.freq_lo)
        self.samples_per_dm = (dt1_f - dt1_i) / self.dt_sample
        
        # idm0, idm1 = min/min dedispersion delay of search, in samples
        self.idm0 = int(sparams.dm_min * self.samples_per_dm)
        self.idm1 = int(sparams.dm_max * self.samples_per_dm) + 1

        # A little trick to get the normalization, by running FDMT on an all-ones array.

        nt0 = round_up_to_power_of_two(self.idm1 + 50)
        a = np.ones((sparams.nfreq, nt0), dtype=np.float32)
        a = bz_fdmt.FDMT(a, self.freq_lo, self.freq_hi, self.idm1, np.float32, Verbose=False)

        self.var = np.copy(a[:,-1])   # 1D array of length ndm

        # Sanity checks, just to sleep better at night.

        assert a.shape == (self.idm1, nt0)
        assert np.all(self.var > 0.0)
        assert np.all(np.abs(self.var - np.array(self.var+0.5, dtype=np.int)) < 1.0e-5)   # variances should all be integers

        for idm in xrange(self.idm1):
            assert np.all(np.abs(a[idm,idm:] - self.var[idm]) < 1.0e-5 * self.var[idm])


    def dedisperse(self, intensity):
        assert intensity.shape == (self.search_params.nfreq, self.search_params.nsamples)

        # Run FDMT!
        a = bz_fdmt.FDMT(intensity, self.freq_lo, self.freq_hi, self.idm1, np.float32, Verbose=False)

        # Now we just want to find the (DM,time) of the most significant pulse.
        max_snr = 0.0
        max_idm = 0.0
        max_it = 0.0

        for idm in xrange(self.idm0, self.idm1):
            it = int(np.argmax(a[idm,idm:])) + idm
            snr = a[idm,it] / self.var[idm]**0.5

            if (idm == 0) or (snr > max_snr):
                max_snr = snr
                max_idm = idm
                max_it = it

        return { 'snr': max_snr, 
                 'dm': max_idm / self.samples_per_dm, 
                 'tfin': it * self.dt_sample }

    
    def jsonize(self):
        return { }


    @staticmethod
    def from_json(j, filename=None):
        return bz_fdmt_dedisperser()

