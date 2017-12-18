import copy
import json

from utils import dispersion_delay


class search_params:
    """
    Represents a set of parameters for the FRB search.
    
       nfreq: number of frequency channels (int)
       nsamples: number of time samples (int)
       dt_sample: length of a time sample in seconds (float)
       freq_lo_MHz: 
    """

    def __init__(self, nfreq, freq_lo_MHz, freq_hi_MHz, nsamples, dt_sample, dm_max,
                 dm_min=0.0, sm_min=0.0, sm_max=0.0, beta_min=0.0, beta_max=0.0,
                 width_min=0.0, width_max=0.0, filename=None):

        self.nfreq = nfreq
        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.nsamples = nsamples
        self.dt_sample = dt_sample
        self.dm_min = dm_min
        self.dm_max = dm_max
        self.sm_min = sm_min
        self.sm_max = sm_max
        self.beta_min = beta_min
        self.beta_max = beta_max
        self.width_min = width_min
        self.width_max = width_max

        # Argument checking follows.
        
        if filename is None:
            filename = 'frb_olympics.search_params constructor'

        assert self.dm_min >= 0.0, filename + ": failed assert: dm_min >= 0.0"
        assert self.sm_min >= 0.0, filename + ": failed assert: sm_min >= 0.0"
        assert self.width_min >= 0.0, filename + ": failed assert: width_min >= 0.0"
        assert self.nsamples > 0, filename + ": failed assert: nsamples > 0"
        assert self.nfreq > 0, filename + ": failed assert: nfreq > 0"

        # The choice of ranges here is intended to guard against accidentally using the wrong units
        # (e.g. GHz instead of MHz, millseconds instead of seconds)
        assert self.freq_lo_MHz >= 100.0, filename + ": failed assert: freq_lo_MHz >= 100.0"
        assert self.dt_sample >= 2.0e-6, filename + ": failed assert: dt_sample >= 2.0e-6"
        assert self.dt_sample <= 0.01, filename + ": failed assert: dt_sample <= 0.01"

        assert self.dm_min <= self.dm_max, filename + ": failed assert: dm_min <= dm_max"
        assert self.sm_min <= self.sm_max, filename + ": failed assert: sm_min <= sm_max"
        assert self.width_min <= self.width_max, filename + ": failed assert: width_min <= width_max"
        assert self.freq_lo_MHz < self.freq_hi_MHz, filename + ": failed assert: freq_lo_MHz < freq_hi_MHz"

        # One last, less trivial consistency test: the total timestream length must be at least 
        # 20% larger than the max in-band dispersion delay.

        timestream_length = self.nsamples * self.dt_sample
        max_dispersion_delay = dispersion_delay(self.dm_max, self.freq_lo_MHz) - dispersion_delay(self.dm_max, self.freq_hi_MHz)

        assert timestream_length > 1.2 * max_dispersion_delay, \
            filename + ': failed assert: timestream_length > 1.2 * max_dispersion_delay'
        

    def jsonize(self):
        return { 'nfreq': self.nfreq,
                 'freq_lo_MHz': self.freq_lo_MHz,
                 'freq_hi_MHz': self.freq_hi_MHz,
                 'nsamples': self.nsamples,
                 'dt_sample': self.dt_sample,
                 'dm_min': self.dm_min,
                 'dm_max': self.dm_max,
                 'sm_min': self.sm_min,
                 'sm_max': self.sm_max,
                 'beta_min': self.beta_min,
                 'beta_max': self.beta_max,
                 'width_min': self.width_min,
                 'width_max': self.width_max }

    
    @staticmethod
    def from_json(j, filename=None):
        if filename is None:
            filename = 'frb_olympics.search_params.from_json()'
        if not isinstance(j, dict):
            raise RuntimeError('%s: expected dict, got %s' % (filename, j.__class__.__name__))

        required_keys = set(['nfreq', 'freq_lo_MHz', 'freq_hi_MHz', 'nsamples', 'dt_sample', 'dm_max'])
        optional_keys = set(['dm_min', 'sm_min', 'sm_max', 'beta_min', 'beta_max', 'width_min', 'width_max'])
        specified_keys = set([ str(x) for x in j.keys() ])
        
        missing_keys = required_keys.difference(specified_keys)
        unrecognized_keys = set(specified_keys).difference(required_keys).difference(optional_keys)

        if len(missing_keys) > 0:
            raise RuntimeError('%s: missing key(s): %s' % (filename, sorted(missing_keys)))
        if len(unrecognized_keys) > 0:
            raise RuntimeError('%s: unrecognized key(s): %s' % (filename, sorted(unrecognized_keys)))

        jj = copy.copy(j)
        jj['filename'] = filename
        
        return search_params(**jj)

    
    @staticmethod
    def from_filename(filename):
        f = open(filename)

        try:
            j = json.load(f)
        except:
            raise RuntimeError("%s: couldn't parse json file" % filename)
        
        return search_params.from_json(j, filename)
