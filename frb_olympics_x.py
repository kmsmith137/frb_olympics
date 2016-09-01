# import rf_pipelines


class search_params:
    def __init__(self, filename):
        self.filename = filename

        # Note: extract_field() calls setattr() to set the attribute by the same name
        kv_pairs = self.parse_file(filename)
        self.extract_field('dm_min', float, kv_pairs)
        self.extract_field('dm_max', float, kv_pairs)
        self.extract_field('sm_min', float, kv_pairs)
        self.extract_field('sm_max', float, kv_pairs)
        self.extract_field('beta_min', float, kv_pairs)
        self.extract_field('beta_max', float, kv_pairs)
        self.extract_field('width_sec_min', float, kv_pairs)
        self.extract_field('width_sec_max', float, kv_pairs)
        self.extract_field('nfreq', int, kv_pairs)
        self.extract_field('freq_lo_MHz', float, kv_pairs)
        self.extract_field('freq_hi_MHz', float, kv_pairs)
        self.extract_field('dt_sec', float, kv_pairs)
        self.extract_field('nsamples', int, kv_pairs)

        if len(kv_pairs) > 0:
            raise RuntimeError("%s: unrecognized parameter(s) in file: %s" % (filename, ', '.join(kv_pairs.keys())))

        assert self.dm_min >= 0.0, filename + ": failed assert: dm_min >= 0.0"
        assert self.sm_min >= 0.0, filename + ": failed assert: sm_min >= 0.0"
        assert self.width_sec_min >= 0.0, filename + ": failed assert: width_sec_min >= 0.0"
        assert self.nsamples > 0, filename + ": failed assert: nsamples > 0"
        assert self.nfreq > 0, filename + ": failed assert: nfreq > 0"

        # The choice of ranges here is intended to guard against accidentally using the wrong units
        # (e.g. GHz instead of MHz, millseconds instead of seconds)
        assert self.freq_lo_MHz >= 100.0, filename + ": failed assert: freq_lo_MHz >= 100.0"
        assert self.dt_sec >= 2.0e-6, filename + ": failed assert: dt_sec >= 2.0e-6"
        assert self.dt_sec <= 0.01, filename + ": failed assert: dt_sec <= 0.01"

        assert self.dm_min <= self.dm_max, filename + ": failed assert: dm_min <= dm_max"
        assert self.sm_min <= self.sm_max, filename + ": failed assert: sm_min <= sm_max"
        assert self.width_sec_min <= self.width_sec_max, filename + ": failed assert: width_sec_min <= width_sec_max"
        assert self.freq_lo_MHz < self.freq_hi_MHz, filename + ": failed assert: freq_lo_MHz < freq_hi_MHz"


    @staticmethod
    def parse_file(filename):
        """Reads key/value pairs from file, and returns a string->string dictionary."""

        ret = { }

        for line in open(filename):
            if (len(line) > 0) and (line[-1] == '\n'):
                line = line[:-1]

            i = line.find('#')
            if i >= 0:
                line = line[:i]

            t = line.split()
            if len(t) == 0:
                continue
            if (len(t) != 3) or (t[1] != '='):
                raise RuntimeError("%s: parse error in line '%s'" % (filename,line))
            if ret.has_key(t[0]):
                raise RuntimeError("%s: duplicate key '%s'" % (filename,t[0]))
            ret[t[0]] = t[2]

        return ret

    
    def extract_field(self, field_name, field_type, kv_pairs):
        if not kv_pairs.has_key(field_name):
            raise RuntimeError("%s: field '%s' not found" % (self.filename, field_name))

        try:
            field_value = field_type(kv_pairs.pop(field_name))
        except:
            raise RuntimeError("%s: parse error in field '%s' (value='%s')" % (self.filename, field_name, kv_pairs[field_name]))

        setattr(self, field_name, field_value)
