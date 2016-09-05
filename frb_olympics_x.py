import sys
import copy
import json
import numpy as np

import rf_pipelines


# Note: cut-and-paste from "simpulse"
def dispersion_delay(dm, freq_MHz):
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);


####################################################################################################


def run_monte_carlo(sp, stream, dedisperser_list, snr=30.):
    """Returns json object (not string)."""

    assert isinstance(sp, search_params)
    assert isinstance(stream, rerunnable_gaussian_noise_stream)
    assert all(isinstance(dedisperser, rf_pipelines.wi_transform) for dedisperser in dedisperser_list)

    # Generate random FRB params
    true_dm = np.random.uniform(sp.dm_min, sp.dm_max)
    true_sm = np.random.uniform(sp.sm_min, sp.sm_max)
    true_beta = np.random.uniform(sp.beta_min, sp.beta_max)
    true_width = np.random.uniform(sp.width_sec_min, sp.width_sec_max)
    
    # Dispersion delays at DM of FRB
    true_dt_initial = dispersion_delay(true_dm, sp.freq_hi_MHz)
    true_dt_final = dispersion_delay(true_dm, sp.freq_lo_MHz)

    timestream_length = sp.nsamples * sp.dt_sec
    tu_min = 0.05*timestream_length - true_dt_initial
    tu_max = 0.95*timestream_length - true_dt_final
    
    assert tu_min < tu_max
    true_tu = np.random.uniform(tu_min, tu_max);
    true_tc = true_tu + (true_dt_initial + true_dt_final) / 2.

    t_frb = rf_pipelines.frb_injector_transform(snr = snr,
                                                undispersed_arrival_time = true_tu,
                                                dm = true_dm,
                                                intrinsic_width = true_width,
                                                sm = true_sm,
                                                spectral_index = true_beta)

    # Start building up output
    json_output = { 'true_snr': snr,
                    'true_dm': true_dm,
                    'true_sm': true_sm,
                    'true_beta': true_beta,
                    'true_width': true_width,
                    'true_tcentral': true_tc,
                    'search_results': [ ] }

    saved_state = stream.get_state()

    for dedisperser in dedisperser_list:
        stream.set_state(saved_state)

        pipeline_json = stream.run([t_frb, dedisperser], return_json=True)
        pipeline_json = json.loads(pipeline_json)   # Temporary kludge (I think)
        pipeline_json = pipeline_json[0]['transforms'][-1]

        if not pipeline_json.has_key('frb_global_max_trigger'):
            raise RuntimeError("internal error: dedisperser transform didn't output 'frb_global_max_trigger' field")
        if not pipeline_json.has_key('frb_global_max_trigger_dm'):
            raise RuntimeError("internal error: dedisperser transform didn't output 'frb_global_max_trigger_dm' field")
        
        recovered_snr = pipeline_json['frb_global_max_trigger']
        recovered_dm = pipeline_json['frb_global_max_trigger_dm']

        recovered_dt_initial = dispersion_delay(recovered_dm, sp.freq_hi_MHz)
        recovered_dt_final = dispersion_delay(recovered_dm, sp.freq_lo_MHz)

        if pipeline_json.has_key('frb_global_max_trigger_tcentral'):
            recovered_tc = pipeline_json['frb_global_max_trigger_tcentral']
        elif pipeline_json.has_key('frb_global_max_trigger_tinitial'):
            recovered_tc = pipeline_json['frb_global_max_trigger_tinitial'] + (recovered_dt_final - recovered_dt_initial) / 2.0
        elif pipeline_json.has_key('frb_global_max_trigger_tfinal'):
            recovered_tc = pipeline_json['frb_global_max_trigger_tfinal'] - (recovered_dt_final - recovered_dt_initial) / 2.0
        elif pipeline_json.has_key('frb_global_max_trigger_final'):
            recovered_tc = pipeline_json['frb_global_max_trigger_tundisp'] - (recovered_dt_initial + recovered_dt_final) / 2.0
        else:
            raise RuntimeError("internal error: dedisperser transform didn't output 'frb_global_max_trigger_t*' field")

        search_results_json = { 'recovered_snr': recovered_snr,
                                'recovered_dm': recovered_dm,
                                'recovered_tcentral' : recovered_tc }

        json_output['search_results'].append(search_results_json)

    return json_output


####################################################################################################


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

        # One last, less trivial consistency test: the total timestream length must be at least 
        # 20% larger than the max in-band dispersion delay.

        timestream_length = self.nsamples * self.dt_sec
        max_dispersion_delay = dispersion_delay(self.dm_max, self.freq_lo_MHz) - dispersion_delay(self.dm_max, self.freq_hi_MHz)

        assert timestream_length > 1.2 * max_dispersion_delay, \
            filename + ': failed assert: timestream_length > 1.2 * max_dispersion_delay'


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


    def make_stream(self, no_noise_flag=False, nt_chunk=None):
        return rerunnable_gaussian_noise_stream(nfreq = self.nfreq, 
                                                nt_tot = self.nsamples, 
                                                freq_lo_MHz = self.freq_lo_MHz, 
                                                freq_hi_MHz = self.freq_hi_MHz, 
                                                dt_sample = self.dt_sec,
                                                no_noise_flag = no_noise_flag,
                                                nt_chunk = nt_chunk)


####################################################################################################


class rerunnable_gaussian_noise_stream(rf_pipelines.py_wi_stream):
    """
    Similar to rf_pipelines.gaussian_noise_stream, but allows the stream to be rerun by saving its state:
    
        s = rerunnable_gaussian_noise_stream(...)
        saved_state = s.get_state()
           # ... run stream ...
        s.set_state(saved_state)
           # ... rerunning stream will give same output ...

    If 'no_noise_flag' is True, then the stream will output zeroes instead of Gaussian random numbers.
    """

    def __init__(self, nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, no_noise_flag=False, state=None, nt_chunk=None):
        if nt_tot <= 0:
            raise RuntimeError('rerunnable_gaussian_noise_stream constructor: nt_tot must be > 0')
        if nt_chunk is None:
            nt_chunk = min(1024, nt_tot)

        rf_pipelines.py_wi_stream.__init__(self, nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, nt_chunk)

        self.no_noise_flag = no_noise_flag
        self.nt_tot = nt_tot
        self.nt_chunk = nt_chunk
        self.set_state(state)


    def stream_body(self, run_state):
        run_state.start_substream(0.0)

        it = 0
        while it < self.nt_tot:
            nt = min(self.nt_tot-it, self.nt_chunk)
            intensity = np.zeros((self.nfreq,nt), dtype=np.float) if self.no_noise_flag else self.state.standard_normal((self.nfreq,nt))
            weights = np.ones((self.nfreq, nt), dtype=np.float)
            run_state.write(intensity, weights)
            it += self.nt_chunk

        run_state.end_substream()
        

    def get_state(self):
        return copy.copy(self.state)


    def set_state(self, state):
        if state is None:
            self.state = np.random.RandomState()
        else:
            assert isinstance(state, np.random.RandomState)
            self.state = copy.copy(state)
