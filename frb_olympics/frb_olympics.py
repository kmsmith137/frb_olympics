import os
import sys
import copy
import json
import itertools
import importlib
import numpy as np

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import rf_pipelines


# Note: cut-and-paste from "simpulse"
def dispersion_delay(dm, freq_MHz):
    """Returns dispersion delay in seconds, for given DM and frequency."""
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);


####################################################################################################


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


####################################################################################################


class dedisperser_base:
    def init_search_params(self, sp):
        """
        The 'sp' argument is an object of type search_params.
        """
        raise RuntimeError("frb_olympics.dedisperser_base.init_search_params() was not overloaded by subclass %s" % self.__class__.__name__)


    def allocate(self):
        pass


    def dedisperse(self, arr):
        """
        The 'arr' argument is a float32 array of shape (nfreq, nsamples).
        """
        raise RuntimeError("frb_olympics.dedisperser_base.dedisperse() was not overloaded by subclass %s" % self.__class__.__name__)

    def deallocate(self):
        pass


    @staticmethod
    def from_json(j, filename=None):
        f = filename if (filename is not None) else 'frb_olympics.dedisperser_base.from_json()' 
        
        try:
            m = importlib.import_module(j['module_name'])
        except ImportError:
            raise ImportError("%s: couldn't import module %s" % (f, j['module_name']))
                              
        c = getattr(m, j['class_name'], None)

        if c is None:
            raise RuntimeError("%s: couldn't find class '%s' in module '%s" % (f, j['class_name'], j['module_name']))
        if not issubclass(c, dedisperser_base):
            raise RuntimeError("%s: expected class %s.%s to be a subclass of frb_olympics.dedisperser_base" % (f, j['module_name'], j['class_name']))
        if c.from_json == dedisperser_base.from_json:
            raise RuntimeError("%s: expected class %s.%s to override dedisperser_base.from_json()" % (f, j['module_name'], j['class_name']))

        f = filename if (filename is not None) else ('%s.%s.from_json()' % (j['module_name'], j['class_name']))
        return c.from_json(j, f)


####################################################################################################


class olympics:
    """
    This is the class used to maintain state for an FRB olympics run.
    See examples/example0_bonsai/run-example0.py for an example of how to use it.
    """

    def __init__(self, sparams, snr=30.0, simulate_noise=True):
        """
        The 'sparams' arg should either be an object of class search_params (see below for definition)
        or a string, assumed to be a search_params filename.

        The 'snr' arg is the signal-to-noise of the FRB, assuming unit variance.
        """

        self.simulate_noise = simulate_noise
        self.snr = snr

        if isinstance(sparams, search_params):
            self.sparams = sparams
        elif isinstance(sparams, basestring):
            # if 'sparams' is a string, then assume it's the filename
            self.sparams = search_params.from_filename(sparams)
        else:
            raise RuntimeError("frb_olympics.olympics constructor: expected 'sparams' argument to be either a string or an object of class frb_olympics.search_params")

        self.stream = rerunnable_gaussian_noise_stream(nfreq = self.sparams.nfreq, 
                                                       nt_tot = self.sparams.nsamples, 
                                                       freq_lo_MHz = self.sparams.freq_lo_MHz, 
                                                       freq_hi_MHz = self.sparams.freq_hi_MHz, 
                                                       dt_sample = self.sparams.dt_sample,
                                                       simulate_noise = simulate_noise)


        # The dedisperser_list is a list of (name, transform) pairs
        # where the 'transform' is an object of class rf_pipelines.pipeline_object.
        self.dedisperser_list = [ ]


    def add_dedisperser(self, transform, name=None):
        """Adds a dedisperser (represented by 'transform', an object of class rf_pipelines.pipeline_object) to the dedisperser_list."""

        assert isinstance(transform, rf_pipelines.pipeline_object)

        if name is None:
            name = transform.name

        self.dedisperser_list.append((name, transform))


    def add_bonsai(self, config_filename, name=None, use_analytic_normalization=True, dm_min=None, dm_max=None):
        """Adds a bonsai_dedisperser to the dedisperser_list."""
 
        t = rf_pipelines.bonsai_dedisperser(config_filename, 
                                            fill_rfi_mask = False, 
                                            track_global_max = True, 
                                            use_analytic_normalization = use_analytic_normalization,
                                            dm_min = dm_min,
                                            dm_max = dm_max)

        self.add_dedisperser(t, name)


    def add_bb_dedisperser(self, dm_tol, dm_t0, name=None, verbosity=1):
        """Adds a bb_dedisperser to the dedisperser_list."""

        t = rf_pipelines.bb_dedisperser(dm_start = self.sparams.dm_min,
                                        dm_end = self.sparams.dm_max,
                                        dm_tol = dm_tol,
                                        dm_t0 = dm_t0,
                                        scrunch = False,
                                        sc_tol = 1.15,
                                        sc_t0 = 0.0,
                                        nt_in = self.sparams.nsamples,
                                        verbosity = verbosity)

        self.add_dedisperser(t, name)


    def add_bz_fdmt(self, name='FDMT'):
        """Adds an FDMT dedisperser to the dedisperser_list."""

        # No free parameters!
        t = rf_pipelines.bz_fdmt_dedisperser(self.sparams.dm_max, self.sparams.nsamples)
        self.add_dedisperser(t, name)


    def run(self, json_filename, nmc, clobber=False):
        """Runs a set of Monte Carlo simulations."""

        if not json_filename.endswith('.json'):
            raise RuntimeError("frb_olympics.olympics.run(): 'json_filename' argument must end in .json")
        
        if len(self.dedisperser_list) == 0:
            raise RuntimeError('frb_olympics.olympics.run(): no dedispersers were defined')

        if nmc <= 0:
            raise RuntimeError('frb_olympics.olympics.run(): expected nmc > 0')

        json_out = { 'search_params': self.sparams.jsonize(),
                     'simulate_noise': self.simulate_noise,
                     'snr': self.snr,
                     'nmc': nmc,
                     'dedisperser_names': [ ],
                     'sims': [ ] }
        
        for (dedisperser_name, transform) in self.dedisperser_list:
            json_out["dedisperser_names"].append(dedisperser_name)

        if not clobber and os.path.exists(json_filename) and (os.stat(json_filename).st_size > 0):
            for i in itertools.count():
                filename2 = '%s.old%d.json' % (json_filename[:-5], i)
                if not os.path.exists(filename2):
                    print >>sys.stderr, 'renaming existing file %s -> %s' % (json_filename, filename2)
                    os.rename(json_filename, filename2)
                    break

        verb = 'truncating' if os.path.exists(json_filename) else 'creating'
        print >>sys.stderr, verb, 'file', json_filename
        f_out = open(json_filename, 'w')

        for imc in xrange(nmc):
            print >>sys.stderr, 'frb_olympics: starting Monte Carlo %d/%d' % (imc, nmc)
            json_sim = self.run_one()
            json_out['sims'].append(json_sim)

        json.dump(json_out, f_out, indent=4)
        print >>f_out   # extra newline (cosmetic)
        print >>sys.stderr, 'wrote', json_filename
        
        # make_snr_plots() is defined later in this file
        make_snr_plots(json_filename, json_out)


    def run_one(self):
        """Runs a single Monte Carlo simulation (helper function called by run())."""
        
        if len(self.dedisperser_list) == 0:
            raise RuntimeError('frb_olympics.olympics.run_one(): no dedispersers were defined')

        # Generate random FRB params
        true_dm = np.random.uniform(self.sparams.dm_min, self.sparams.dm_max)
        true_sm = np.random.uniform(self.sparams.sm_min, self.sparams.sm_max)
        true_beta = np.random.uniform(self.sparams.beta_min, self.sparams.beta_max)
        true_width = np.random.uniform(self.sparams.width_min, self.sparams.width_max)
    
        # Dispersion delays at DM of FRB
        true_dt_initial = dispersion_delay(true_dm, self.sparams.freq_hi_MHz)
        true_dt_final = dispersion_delay(true_dm, self.sparams.freq_lo_MHz)

        # Min/max allowed _undispersed_ arrival time
        timestream_length = self.sparams.nsamples * self.sparams.dt_sample
        tu_min = 0.05*timestream_length - true_dt_initial
        tu_max = 0.95*timestream_length - true_dt_final
        
        assert tu_min < tu_max
        true_tu = np.random.uniform(tu_min, tu_max);

        # Convert undispersed arrival time to central arrival time.
        true_tc = true_tu + (true_dt_initial + true_dt_final) / 2.
        
        t_frb = rf_pipelines.frb_injector_transform(snr = self.snr,
                                                    undispersed_arrival_time = true_tu,
                                                    dm = true_dm,
                                                    variance = 1.0,
                                                    intrinsic_width = true_width,
                                                    sm = true_sm,
                                                    spectral_index = true_beta)

        # Start building up output
        json_output = { 'true_snr': self.snr,
                        'true_dm': true_dm,
                        'true_sm': true_sm,
                        'true_beta': true_beta,
                        'true_width': true_width,
                        'true_tcentral': true_tc,
                        'search_results': [ ] }

        # We save the RNG state and restore it below, so that each dedisperser "sees" the same
        # noise realization.  At the end of run_one(), the RNG state has been advanced, so that
        # subsequent calls to run_one() will produce a different noise realization.
        saved_state = self.stream.get_state()

        for (name, dedisperser) in self.dedisperser_list:
            print >>sys.stderr, 'frb_olympics: running dedisperser', name
            self.stream.set_state(saved_state)

            p = rf_pipelines.pipeline([self.stream, t_frb, dedisperser])
            pipeline_json = p.run(outdir=None)
            p.unbind()
            
            # We're only interested in the part of the json output from the last transform (the dedisperser).
            pipeline_json = pipeline_json['pipeline'][-1]

            if not pipeline_json.has_key('frb_global_max_trigger'):
                raise RuntimeError("internal error: dedisperser transform didn't output 'frb_global_max_trigger' field")
            if not pipeline_json.has_key('frb_global_max_trigger_dm'):
                raise RuntimeError("internal error: dedisperser transform didn't output 'frb_global_max_trigger_dm' field")
        
            recovered_snr = pipeline_json['frb_global_max_trigger']
            recovered_dm = pipeline_json['frb_global_max_trigger_dm']
            
            recovered_dt_initial = dispersion_delay(recovered_dm, self.sparams.freq_hi_MHz)
            recovered_dt_final = dispersion_delay(recovered_dm, self.sparams.freq_lo_MHz)

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

            print >>sys.stderr, 'frb_olympics: dedisperser=%s, recovered snr=%s' % (name, recovered_snr)

            search_results_json = { 'recovered_snr': recovered_snr,
                                    'recovered_dm': recovered_dm,
                                    'recovered_tcentral': recovered_tc }

            json_output['search_results'].append(search_results_json)
            
        return json_output


####################################################################################################


def make_snr_plot(plot_filename, xvec, snr_arr, xmin, xmax, xlabel, dedisperser_names):
    """Helper function called by make_snr_plots()."""

    color_list = [ 'b', 'r', 'm', 'g', 'k', 'y', 'c' ]

    xvec = np.array(xvec)
    nds = len(dedisperser_names)
    nmc = len(xvec)

    assert xvec.ndim == 1
    assert nds <= len(color_list)
    assert snr_arr.shape == (nmc, nds)
    assert np.all(snr_arr >= 0.0)

    slist = [ ]
    for ids in xrange(nds):
        slist.append(plt.scatter(xvec, snr_arr[:,ids], s=5, color=color_list[ids], marker='o'))

    plt.xlim(xmin, xmax)
    plt.ylim(0.0, max(np.max(snr_arr),1.1))
    plt.axhline(y=1, ls=':', color='k')
    plt.xlabel(xlabel)
    plt.ylabel('Optimality')
    plt.legend(slist, dedisperser_names, scatterpoints=1, loc='lower left')

    print >>sys.stderr, 'writing', plot_filename
    plt.savefig(plot_filename)
    plt.clf()
    

def make_snr_plots(json_filename, json_obj=None):
    """
    Makes plots of snr versus dm.  Also makes plots of snr versus sm/beta/width (if these parameters span a nonzero range).
    This routine can be called to generate plots directly from a json file.
    """

    if not json_filename.endswith('.json'):
        raise RuntimeError("frb_olympics.make_snr_plots(): 'json_filename' argument must end in .json")

    if json_obj is None:
        print >>sys.stderr, 'reading', json_filename
        json_obj = json.load(open(json_filename))

    dedisperser_names = json_obj["dedisperser_names"];
    nmc = json_obj["nmc"]
    nds = len(dedisperser_names)

    snr_arr = np.zeros((nmc,nds))
    for imc in xrange(nmc):
        true_snr = json_obj['sims'][imc]['true_snr']
        for ids in xrange(nds):
            recovered_snr = json_obj['sims'][imc]['search_results'][ids]['recovered_snr']
            snr_arr[imc,ids] = recovered_snr / true_snr

    make_snr_plot(plot_filename = ('%s_snr_vs_dm.pdf' % json_filename[:-5]),
                  xvec = [ s['true_dm'] for s in json_obj['sims'] ],
                  snr_arr = snr_arr,
                  xmin = json_obj['search_params']['dm_min'],
                  xmax = json_obj['search_params']['dm_max'],
                  xlabel = 'DM',
                  dedisperser_names = dedisperser_names)

    if json_obj['search_params']['sm_min'] < json_obj['search_params']['sm_max']:
        make_snr_plot(plot_filename = ('%s_snr_vs_sm.pdf' % json_filename[:-5]),
                      xvec = [ s['true_sm'] for s in json_obj['sims'] ],
                      snr_arr = snr_arr,
                      xmin = json_obj['search_params']['sm_min'],
                      xmax = json_obj['search_params']['sm_max'],
                      xlabel = 'SM',
                      dedisperser_names = dedisperser_names)

    if json_obj['search_params']['beta_min'] < json_obj['search_params']['beta_max']:
        make_snr_plot(plot_filename = ('%s_snr_vs_beta.pdf' % json_filename[:-5]),
                      xvec = [ s['true_beta'] for s in json_obj['sims'] ],
                      snr_arr = snr_arr,
                      xmin = json_obj['search_params']['beta_min'],
                      xmax = json_obj['search_params']['beta_max'],
                      xlabel = 'spectral index',
                      dedisperser_names = dedisperser_names)

    if json_obj['search_params']['width_min'] < json_obj['search_params']['width_max']:
        make_snr_plot(plot_filename = ('%s_snr_vs_width.pdf' % json_filename[:-5]),
                      xvec = [ s['true_width'] for s in json_obj['sims'] ],
                      snr_arr = snr_arr,
                      xmin = json_obj['search_params']['width_min'],
                      xmax = json_obj['search_params']['width_max'],
                      xlabel = 'intrinsic width [sec]',
                      dedisperser_names = dedisperser_names)


####################################################################################################


class rerunnable_gaussian_noise_stream(rf_pipelines.wi_stream):
    """
    Similar to rf_pipelines.gaussian_noise_stream, but allows the stream to be rerun by saving its state:
    
        s = rerunnable_gaussian_noise_stream(...)
        saved_state = s.get_state()
           # ... run stream ...
        s.set_state(saved_state)
           # ... rerunning stream will give same output ...

    If 'no_noise_flag' is True, then the stream will output zeroes instead of Gaussian random numbers.
    """

    def __init__(self, nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, simulate_noise=True, state=None, nt_chunk=None):
        rf_pipelines.wi_stream.__init__(self, 'rereunnable_gaussian_noise_stream')

        if nt_tot <= 0:
            raise RuntimeError('rerunnable_gaussian_noise_stream constructor: nt_tot must be > 0')
        if nt_chunk is None:
            nt_chunk = min(1024, nt_tot)

        # These are members of the rf_pipelines.wi_stream base class.
        self.nfreq = nfreq
        self.nt_chunk = nt_chunk

        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.dt_sample = dt_sample
        self.simulate_noise = simulate_noise
        self.nt_tot = nt_tot
        self.set_state(state)


    def _bind_stream(self, json_attrs):
        json_attrs['freq_lo_MHz'] = self.freq_lo_MHz
        json_attrs['freq_hi_MHz'] = self.freq_hi_MHz
        json_attrs['dt_sample'] = self.dt_sample


    def _fill_chunk(self, intensity, weights, pos):
        assert intensity.shape == weights.shape == (self.nfreq, self.nt_chunk)

        if self.simulate_noise:
            intensity[:,:] = self.state.standard_normal((self.nfreq, self.nt_chunk))
        else:
            intensity[:,:] = np.zeros((self.nfreq, self.nt_chunk), dtype=np.float32)

        weights[:,:] = np.ones((self.nfreq, self.nt_chunk), dtype=np.float32)

        return (pos + self.nt_chunk) < self.nt_tot
        

    def get_state(self):
        """Returns the current RNG state."""
        return copy.copy(self.state)


    def set_state(self, state):
        """Restores the RNG state to a previous value, obtained by calling get_state()."""

        if state is None:
            self.state = np.random.RandomState()
        else:
            assert isinstance(state, np.random.RandomState)
            self.state = copy.copy(state)
