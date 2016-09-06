import os
import sys
import copy
import json
import itertools
import numpy as np

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import rf_pipelines


# Note: cut-and-paste from "simpulse"
def dispersion_delay(dm, freq_MHz):
    return 4.148806e3 * dm / (freq_MHz * freq_MHz);


####################################################################################################
#
# MPI configuration


is_mpi = True

try:
    import mpi4py
except ImportError:
    is_mpi = False


if is_mpi:
    import mpi4py.MPI
    mpi_rank = mpi4py.MPI.COMM_WORLD.rank
    mpi_size = mpi4py.MPI.COMM_WORLD.size
    mpi_bcast = mpi4py.MPI.COMM_WORLD.bcast
    mpi_gather = lambda x: mpi4py.MPI.COMM_WORLD.gather(x, root=0)
    mpi_barrier = mpi4py.MPI.COMM_WORLD.Barrier
    mpi_rank_within_node = 0    # determined below
    mpi_tasks_per_node = 0      # determined below

    import socket
    my_hostname = socket.gethostname()
    all_hostnames = mpi4py.MPI.COMM_WORLD.allgather(my_hostname)
    mult = { }

    for (i,h) in enumerate(all_hostnames):
        mult[h] = mult.get(h,0) + 1
        if h == my_hostname:
            mpi_tasks_per_node += 1
            if i < mpi_rank:
                mpi_rank_within_node += 1

    # check
    for (k,v) in mult.iteritems():
        if v != mpi_tasks_per_node:
            raise RuntimeError('expected all hostnames to have equal numbers of MPI tasks')
        
else:
    # Hmm, mpi4py failed to load... is this because we're a serial job, or for some other reason?
    suspicious_vars = [ 'OMPI_COMM_WORLD_SIZE',   # openmpi-1.3
                        'SLURM_NPROCS',           # slurm v14
                        'PMI_SIZE' ]              # mpich 2

    for x in suspicious_vars:
        if os.environ.has_key(x):
            print >>sys.stderr, 'Environment variable %s is set, but mpi4py failed to import.'
            print >>sys.stderr, 'This probably means that something is wrong with the mpi4py installation.'
            sys.exit(1)

    # OK, satisfied that we're a serial job!
    mpi_rank = 0
    mpi_size = 1
    mpi_bcast = lambda x: x
    mpi_gather = lambda x: [x]
    mpi_rank_within_node = 0
    mpi_tasks_per_node = 1
    mpi_barrier = lambda: None


def init_mpi_log_files(stem):
    if mpi_size == 1:
        return

    # redirect stdout, stderr to log files which are unique to this MPI rank

    log_filename = '%s.stdout%d' % (stem, mpi_rank)
    sys.stderr.write('redirecting stdout -> %s\n' % log_filename)
    f_log = open(log_filename, 'w')
    os.dup2(f_log.fileno(), sys.stdout.fileno())

    log_filename = '%s.stderr%d' % (stem, mpi_rank)
    sys.stderr.write('redirecting stderr -> %s\n' % log_filename)
    f_log = open(log_filename, 'w')
    os.dup2(f_log.fileno(), sys.stderr.fileno())


####################################################################################################


class olympics:
    def __init__(self, sparams, snr=30.0, simulate_noise=True):
        self.simulate_noise = simulate_noise
        self.snr = snr

        if isinstance(sparams, search_params):
            self.sparams = sparams
        elif isinstance(sparams, basestring):
            # if 'sparams' is a string, then assume it's the filename
            self.sparams = search_params(sparams)
        else:
            raise RuntimeError("frb_olympics.olympics constructor: expected 'sparams' argument to be either a string or an object of class frb_olympics.search_params")

        self.stream = rerunnable_gaussian_noise_stream(nfreq = self.sparams.nfreq, 
                                                       nt_tot = self.sparams.nsamples, 
                                                       freq_lo_MHz = self.sparams.freq_lo_MHz, 
                                                       freq_hi_MHz = self.sparams.freq_hi_MHz, 
                                                       dt_sample = self.sparams.dt_sec,
                                                       simulate_noise = simulate_noise)


        # The dedisperser_list is a list of (name, transform) pairs
        self.dedisperser_list = [ ]


    def add_bonsai(self, config_hdf5_filename, name=None):
        transform = rf_pipelines.bonsai_dedisperser(config_hdf5_filename)

        if name is None:
            name = transform.name

        self.dedisperser_list.append((name, transform))


    def run(self, json_filename, nmc, clobber=False, mpi_log_files=True):
        if not json_filename.endswith('.json'):
            raise RuntimeError("frb_olympics.olympics.run(): 'json_filename' argument must end in .json")
        
        if len(self.dedisperser_list) == 0:
            raise RuntimeError('frb_olympics.olympics.run(): no dedispersers were defined')

        if nmc <= 0:
            raise RuntimeError('frb_olympics.olympics.run(): expected nmc > 0')
        
        if (mpi_size > 1) and mpi_log_files:
            init_mpi_log_files(stem = json_filename[:-5])

        if mpi_rank == 0:
            json_out = { 'simulate_noise': self.simulate_noise,
                         'snr': self.snr,
                         'nmc': nmc,
                         'dedisperser_names': [ ],
                         'search_params': { },
                         'sims': [ ] }
        
            for (field_name, field_type) in self.sparams.all_fields:
                json_out["search_params"][field_name] = getattr(self.sparams, field_name)
            
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

        json_sims = [ ]

        for imc in xrange(mpi_rank, nmc, mpi_size):
            print >>sys.stderr, 'frb_olympics: starting Monte Carlo %d/%d' % (imc, nmc)
            json_sim = self.run_one()
            json_sims.append(json_sim)

        json_sims = mpi_gather(json_sims)

        if mpi_rank == 0:
            for s in json_sims:
                json_out['sims'] += s

            json.dump(json_out, f_out, indent=4)
            print >>f_out   # extra newline (cosmetic)
            print >>sys.stderr, 'wrote', json_filename

            # make_snr_plots() is defined later in this file
            make_snr_plots(json_filename, json_out)


    def run_one(self):
        """Returns json object (not string representation).  Uses current state of stream"""
        
        if len(self.dedisperser_list) == 0:
            raise RuntimeError('frb_olympics.olympics.run_one(): no dedispersers were defined')

        # Generate random FRB params
        true_dm = np.random.uniform(self.sparams.dm_min, self.sparams.dm_max)
        true_sm = np.random.uniform(self.sparams.sm_min, self.sparams.sm_max)
        true_beta = np.random.uniform(self.sparams.beta_min, self.sparams.beta_max)
        true_width = np.random.uniform(self.sparams.width_sec_min, self.sparams.width_sec_max)
    
        # Dispersion delays at DM of FRB
        true_dt_initial = dispersion_delay(true_dm, self.sparams.freq_hi_MHz)
        true_dt_final = dispersion_delay(true_dm, self.sparams.freq_lo_MHz)

        timestream_length = self.sparams.nsamples * self.sparams.dt_sec
        tu_min = 0.05*timestream_length - true_dt_initial
        tu_max = 0.95*timestream_length - true_dt_final
        
        assert tu_min < tu_max
        true_tu = np.random.uniform(tu_min, tu_max);
        true_tc = true_tu + (true_dt_initial + true_dt_final) / 2.
        
        t_frb = rf_pipelines.frb_injector_transform(snr = self.snr,
                                                    undispersed_arrival_time = true_tu,
                                                    dm = true_dm,
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
        
        saved_state = self.stream.get_state()

        for (name, dedisperser) in self.dedisperser_list:
            self.stream.set_state(saved_state)
            
            pipeline_json = self.stream.run([t_frb, dedisperser], outdir=None, return_json=True)
            
            # A kludge: eventually, the run() return value will be a json object, but for now it returns
            # the string representation, which can be converted to a json object by calling json.loads().
            pipeline_json = json.loads(pipeline_json)
            
            # We're only interested in the part of the json output from the last transform (the dedisperser).
            pipeline_json = pipeline_json[0]['transforms'][-1]
            
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
    plt.ylabel('Recovered SNR')
    plt.legend(slist, dedisperser_names, scatterpoints=1, loc='lower left')

    print >>sys.stderr, 'writing', plot_filename
    plt.savefig(plot_filename)
    plt.clf()
    

def make_snr_plots(json_filename, json_obj=None):
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

    if json_obj['search_params']['width_sec_min'] < json_obj['search_params']['width_sec_max']:
        make_snr_plot(plot_filename = ('%s_snr_vs_width.pdf' % json_filename[:-5]),
                      xvec = [ s['true_width'] for s in json_obj['sims'] ],
                      snr_arr = snr_arr,
                      xmin = json_obj['search_params']['width_sec_min'],
                      xmax = json_obj['search_params']['width_sec_max'],
                      xlabel = 'intrinsic width [sec]',
                      dedisperser_names = dedisperser_names)


####################################################################################################



def parse_kv_file(filename):
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


class search_params:
    # a list of (field_name, field_type) pairs
    all_fields = [ ('dm_min', float),
                   ('dm_max', float),
                   ('sm_min', float),
                   ('sm_max', float),
                   ('beta_min', float),
                   ('beta_max', float),
                   ('width_sec_min', float),
                   ('width_sec_max', float),
                   ('nfreq', int),
                   ('freq_lo_MHz', float),
                   ('freq_hi_MHz', float),
                   ('dt_sec', float),
                   ('nsamples', int) ]


    def __init__(self, filename):
        kv_pairs = parse_kv_file(filename)

        for (field_name, field_type) in self.all_fields:
            try:
                field_value = kv_pairs.pop(field_name)
                field_value = field_type(field_value)
            except KeyError:
                raise RuntimeError("%s: field '%s' not found" % (filename, field_name))
            except ValueError:
                raise RuntimeError("%s: parse error in field '%s' (value='%s')" % (filename, field_name, field_value))

            setattr(self, field_name, field_value)

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

    def __init__(self, nfreq, nt_tot, freq_lo_MHz, freq_hi_MHz, dt_sample, simulate_noise=True, state=None, nt_chunk=None):
        if nt_tot <= 0:
            raise RuntimeError('rerunnable_gaussian_noise_stream constructor: nt_tot must be > 0')
        if nt_chunk is None:
            nt_chunk = min(1024, nt_tot)

        rf_pipelines.py_wi_stream.__init__(self, nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, nt_chunk)

        self.simulate_noise = simulate_noise
        self.nt_tot = nt_tot
        self.nt_chunk = nt_chunk
        self.set_state(state)


    def stream_body(self, run_state):
        run_state.start_substream(0.0)

        it = 0
        while it < self.nt_tot:
            nt = min(self.nt_tot-it, self.nt_chunk)
            intensity = self.state.standard_normal((self.nfreq,nt)) if self.simulate_noise else np.zeros((self.nfreq,nt), dtype=np.float)
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
