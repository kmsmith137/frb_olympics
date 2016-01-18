"""
Here are some notes on the frb_olympics internals (see also the README file)

One important data structure is frb_search_params, which is shared
between C++ and Python, and defines various parameters of the search.

The frb_search_params constructor syntax is
   p = frb_search_params('example_search_params.txt')    # construct from file
   p = frb_search_params(p2)           # make copy from existing search_params

An frb_search_params contains the following members

   # These fields define the parameter space to be simulated or searched
   p.dm_min, p.dm_max         # dispersion measure (DM) range
   p.sm_min, p.sm_max         # scattering measure (SM) range, not really used yet
   p.beta_min, p.beta_max     # spectral index range, not really used yet
   p.width_min, p.width_max   # intrinsic width range

   # These fields define the frequency channels of the instrument, and its time
   # resolution.  Channels are assumed equally spaced in frequency nu.
   p.nchan               # number of frequency channels
   p.band_lo_freq_MHz    # e.g. 400.0 for full CHIME
   p.band_hi_freq_MHz    # e.g. 800.0 for full CHIME
   p.dt_sample           # e.g. 1.0e-3 for full CHIME (all times are in seconds)

   # These fields define the length of the timestream to be searched.
   # For a non-incremental search, we have nsamples_tot == nsamples_per_chunk and nchunks == 1.
   # For an incremental search, we have nsamples_tot == nsamples_per_chunk * nchunks.
   p.nsamples_tot
   p.nsamples_per_chunk
   p.nchunks

For the full declaration of frb_search_params, see frb_olympics.hpp (C++) or 
frb_olympics_c.pyx (Python).  The python declaration is sort of hard to read 
due to all the cython boilerplate!

I usually initialize the frb_search_params by writing a .txt file (best illustrated 
by example; see example_search_params.txt) and then doing:
   p = frb_olympics.frb_search_params('example_search_params.txt')


The frb_olympics uses a plugin architecture in which new search algorithms can
be included as long as they implement an appropriate API.  For search algorithms
written in Python, the API is as follows.  Each search algorithm should be a class
which implements the following methods:

  __init__(self, ...):
      should initialize self.name, in addition to any fields needed internally

  search_init(self, search_params):
      Here, 'search_params' is an object of class frb_search_params, see above.
      This method is called as soon as the search params are determined.  It will
      only be called once, even if multiple timestreams are searched.

      It should initialize the following fields, in addition to any fields
      needed internally:

         # Memory used while searching, in GB
         self.search_gb

         # Normally, the search algorithm will just compute the trigger statistic (a scalar)
         # without actually returning the output of its DM transform (a 2D array indexed by 
         # DM and arrival time).  
         #
         # However, for debugging purposes we sometimes want to inspect the output of the transform.
         # The algorithm should set (self.debug_buffer_ndm, self.debug_buffer_nt) to the shape
         # of the transform output array.
         #
         self.debug_buffer_ndm
         self.debug_buffer_nt

   search_start(self, mpi_rank_within_node):
       Called once per searched timestream (unlike search_init(), which is only called once, period).
       Does per-timestream initializations, such as allocating buffers.
       The 'mpi_rank_within_node' argument is sometimes useful for pinning threads to cores.

   search_chunk(self, chunk, ichunk, debug_buffer=None):
       This is the main routine which does the actual search.  Its job is to set
       self.search_result to the "trigger" statistic T from the FRB olympics memo.  

       In a non-incremental search, search_chunk() will be called once with ichunk=0.

       In an incremental search, search_chunk() will be called multiple times (with 
       ichunk=0,..,nchunks-1) and search_chunk() is responsible for saving state between
       calls and setting self.search_result at the end.

       Normally, search_chunk() just computes the trigger statistic (a scalar) without
       actually returning the output of the DM transform (a 2D array indexed by DM and
       arrival time).  However, if debug_buffer is not None, it will be an array of shape
       (self.debug_buffer_ndm, self.debug_buffer_nt) which holds the output of the DM
       transform.  I usually use this for visual inspection with frb-dump.py, but it
       could be used for other things.  

       "Masked" entries in the debug_buffer which are not actually searched
       (because their DM or arrival time is out of range) are indicated by setting
       them to -1.0e30.

  search_end(self):
       This is called once per searched timestream, after the calls to search_chunk().
       Usually it just deallocates buffers, which is important so that we don't use
       too much memory when multiple searches a run in parallel with frb-compare.py.
       

Summarizing, the structure of frb-compare.py is approximately this:

   for each search algorithm a:
      call a.__init__()
      call a.search_init(search_params)

   for each timestream:
      for each algorithm a:
         call a.search_start()
         for each chunk in an incremental search:
             call a.search_chunk()
         call a.search_end()

For a well-commented example of implementing the search algorithm API in Python,
see frb_fdmt.py (just a wrapper around Alex Josephy's cpuFDMT.py, but shows how 
to implement the API).

See the README file for some frb-compare.py example runs.  Another useful utility
for visual sanity checking is frb-dump.py, see the bottom of the README file.
"""

import os
import sys
import numpy as np

# cython imports, including the search algorithms 
#   { simple_direct, sloth, bonsai }
from frb_olympics_c import frb_rng, frb_pulse, frb_search_params, \
    simple_direct, sloth, bonsai, sloth_sm_subsearch

# Alex Josephy FDMT code
from frb_fdmt import fdmt

# the downsample search "algorithm", which just wraps
# another search algorithm and runs it at lower time sampling
from frb_downsample import downsample

# the rechunk search "algorithm" wraps another search algorithm and runs it
# with a different chunk size (useful for debugging incremental search)
from frb_rechunk import rechunk

# another wrapper algorithm; this one combines multiple algorithms
# (e.g. trees with different dm_max and downsampling) and returns the max trigger
from frb_combiner import combiner

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

    # redirect stdout to log file which is unique to this MPI rank
    log_filename = '%s_log%d' % (stem, mpi_rank)
    print 'redirecting stdout -> %s' % log_filename

    f_log = open(log_filename, 'w')
    os.dup2(f_log.fileno(), sys.stdout.fileno())


def imp(filename):
    """
    A utility routine which imports module with given filename, returning the module object.

    In the MPI case, this must be called MPI-collectively, and is written to
    avoid a race condition where one core is writing the .pyc file while another
    core concurrently reads the incomplete .pyc file.
    """

    import imp

    if mpi_rank > 0:
        e = mpi_bcast(None)
        if e is not None:
            raise e

    try:
        module_name = os.path.basename(filename)
        module_name = module_name[:module_name.find('.')]
        ret = imp.load_source(module_name, filename)
    except Exception as e:
        if mpi_rank == 0:
            mpi_bcast(e)
        raise

    if mpi_rank == 0:
        mpi_bcast(None)

    return ret


####################################################################################################
#
# Algorithm registry


algo_list = [ ]
memhack_list = [ ]
ini_flag = False


def add_algo(algo, memhack=1):
    assert memhack >= 1
    algo_list.append(algo)
    memhack_list.append(memhack)


def init_algorithms(search_params):
    global ini_flag

    if ini_flag:
        raise RuntimeError('double call to frb_olympics.init_algorithms()')

    assert len(algo_list) > 0

    for algo in algo_list:
        algo.search_init(search_params)

        required_fields = [ 'name', 'search_params', 'debug_buffer_ndm', 
                            'debug_buffer_nt', 'search_gb', 'search_init', 
                            'search_start', 'search_chunk', 'search_end' ]

        missing_fields = [ f for f in required_fields if not hasattr(algo,f) ]

        if len(missing_fields) > 0:
            raise RuntimeError('algorithm object is missing the following required fields: %s' % missing_fields)

        assert isinstance(algo.name, basestring)
        assert len(algo.name) > 0
        assert algo.debug_buffer_ndm > 0
        assert algo.debug_buffer_nt > 0
        assert algo.search_gb >= 0.0

    ini_flag = True


def enumerate_algorithms_with_memhack():
    """Generates (ialgo,algo) pairs.  Note that this calls search_start() and search_end() "under the hood"."""

    for (ialgo,(algo,memhack)) in enumerate(zip(algo_list,memhack_list)):
        assert memhack > 0
        assert mpi_tasks_per_node % memhack == 0

        nbarriers1 = mpi_rank_within_node % memhack
        nbarriers2 = memhack - nbarriers1

        for i in xrange(nbarriers1):
            mpi_barrier()

        algo.search_start(mpi_rank_within_node // memhack)

        yield (ialgo, algo)

        algo.search_end()
        
        for i in xrange(nbarriers2):
            mpi_barrier()


####################################################################################################


class comparison_outputs:
    def __init__(self, noise_data, pulse_data, search_params=None):
        assert noise_data.ndim == 2
        assert pulse_data.ndim == 2
        assert pulse_data.shape[1] == noise_data.shape[1] + 5

        self.npulse = pulse_data.shape[0]
        self.nnoise = noise_data.shape[0]
        self.nalgo = noise_data.shape[1]

        if self.nnoise >= 2:
            self.noise_mean = np.mean(noise_data, axis=0)
            self.noise_stddev = np.std(noise_data, axis=0)

        self.pulse_arrival_time = pulse_data[:,0]
        self.pulse_intrinsic_width = pulse_data[:,1]
        self.pulse_dm = pulse_data[:,2]
        self.pulse_sm = pulse_data[:,3]
        self.pulse_beta = pulse_data[:,4]
        self.pulse_data = pulse_data[:,5:]
        self.noise_data = noise_data

        if self.nnoise >= 2 and self.npulse >= 1:
            # "Sigma" statistic (to be replaced later by something better)
            self.pulse_sigma = (self.pulse_data - self.noise_mean[np.newaxis,:]) / self.noise_stddev[np.newaxis,:]

        self.search_params = search_params


    def plot_histogram(self, ialgo, filename, xmax=None):
        import matplotlib.pyplot as plt

        assert self.nnoise > 0 and self.npulse > 0
        assert 0 <= ialgo < self.nalgo

        plt.hist(self.noise_data[:,ialgo], label='noise-only sims')
        plt.hist(self.pulse_data[:,ialgo], label='pulse-only sims')

        plt.xlim(xmin=0)
        if xmax is not None:
            plt.xlim(xmax=xmax)
        
        plt.xlabel('$T$')
        plt.ylabel('Counts')
        plt.legend()

        plt.savefig(filename)
        plt.clf()
        print 'wrote', filename


    def plot_sigma(self, stem, ialgo_list=None, color_list=None, marker_list=None, label_list=None, legloc='lower left', **kwds):
        """
        Writes two files ${stem}_sigma.pdf and ${stem}_sigmap.pdf

        The **kwds can be any of: dm_min, dm_max, sm_min, sm_max, beta_min, beta_max

        These influence the plot xranges, which are given by (in order of priority):
           - kwds to this routine
           - values in the search_params specified at construction, if not None
           - matplotlib defaults
        """

        import matplotlib.pyplot as plt

        assert self.npulse > 0

        if ialgo_list is None:
            ialgo_list = range(self.nalgo)

        if color_list is None:
            color_list = [ 'b', 'r', 'm', 'g', 'k', 'y', 'c' ]

        if marker_list is None:
            marker_list = [ 'o' for i in xrange(len(ialgo_list)) ]

        if label_list is None:
            assert len(algo_list) == self.nalgo
            label_list = [ algo.name for algo in algo_list ]

        # replace underscores with hyphens so latex doesn't get confused
        label_list = [ l.replace('_','-') for l in label_list ]

        assert len(color_list) >= len(ialgo_list)
        assert len(marker_list) >= len(ialgo_list)
        assert len(label_list) == len(ialgo_list)

        # (xdata, xlabel, xstem, xmin_str, xmax_str) triples
        x_todo = [ (self.pulse_dm, 'DM', '', 'dm_min', 'dm_max') ]
        if (np.max(self.pulse_sm) - np.min(self.pulse_sm)) > 1.0e-3:
            x_todo.append((self.pulse_sm, 'SM', '_sm', 'sm_min', 'sm_max'))
        if (np.max(self.pulse_beta) - np.min(self.pulse_beta)) > 1.0e-3:
            x_todo.append((self.pulse_beta, 'spectral index', '_beta', 'beta_min', 'beta_max'))
        if (np.max(self.pulse_intrinsic_width) - np.min(self.pulse_intrinsic_width)) > 1.0e-5:
            x_todo.append((self.pulse_intrinsic_width, 'intrinsic_width', '_width', 'width_min', 'width_max'))

        # (ydata, ylabel, ystem) triples
        y_todo = [ (self.pulse_data / 30.0, r"$\sigma'$", '_sigmap') ]
        if hasattr(self, 'pulse_sigma'):
            y_todo.append((self.pulse_sigma, r"$\sigma$", '_sigma'))

                   
        for (xdata, xlabel, xstem, xmin_str, xmax_str) in x_todo:
            for (ydata, ylabel, ystem) in y_todo:
                filename = ('%s%s%s.pdf' % (stem, xstem, ystem))

                slist = [ ]
                for (i,ialgo) in enumerate(ialgo_list):
                    slist.append(plt.scatter(xdata, ydata[:,ialgo], s=5, color=color_list[i], marker=marker_list[i]))

                if kwds.has_key(xmin_str):
                    plt.xlim(xmin = kwds[xmin_str])
                elif self.search_params is not None:
                    plt.xlim(xmin = getattr(self.search_params, xmin_str))

                if kwds.has_key(xmax_str):
                    plt.xlim(xmax = kwds[xmax_str])
                elif self.search_params is not None:
                    plt.xlim(xmax = getattr(self.search_params, xmax_str))

                plt.ylim(ymin=0)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                plt.legend(slist, label_list, scatterpoints=1, loc=legloc)

                plt.savefig(filename)
                plt.clf()
                print 'wrote', filename
