import os
import sys
import numpy as np

# cython imports, including the search algorithms 
#   { simple_direct, simple_tree, sloth, bonsai}
from frb_olympics_c import frb_rng, frb_pulse, frb_search_params, \
    simple_direct, simple_tree, sloth, bonsai

# the downsample search "algorithm", which just wraps
# another search algorithm and runs it at lower time sampling
from frb_downsample import downsample


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
    except e:
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


def enumerate_algorithms_with_memhack(bracket_search=True):
    """Generates (ialgo,algo) pairs."""

    for (ialgo,(algo,memhack)) in enumerate(zip(algo_list,memhack_list)):
        assert memhack > 0
        assert mpi_tasks_per_node % memhack == 0

        nbarriers1 = mpi_rank_within_node % memhack
        nbarriers2 = memhack - nbarriers1 - 1

        for i in xrange(nbarriers1):
            mpi_barrier()

        if bracket_search:
            algo.search_start()

        yield (ialgo, algo)

        if bracket_search:
            algo.search_end()
        
        for i in xrange(nbarriers2):
            mpi_barrier()


####################################################################################################


class comparison_outputs:
    def __init__(self, noise_data, pulse_data):
        assert noise_data.ndim == 2
        assert pulse_data.ndim == 2
        assert pulse_data.shape[1] == noise_data.shape[1] + 5
        assert noise_data.shape[0] >= 2
        assert pulse_data.shape[1] >= 1

        self.noise_data = noise_data
        self.nnoise = self.noise_data.shape[0]
        self.nalgo = self.noise_data.shape[1]
        self.noise_mean = np.mean(self.noise_data, axis=0)
        self.noise_stddev = np.std(self.noise_data, axis=0)

        self.npulse = pulse_data.shape[0]
        self.pulse_arrival_time = pulse_data[:,0]
        self.pulse_intrinsic_width = pulse_data[:,1]
        self.pulse_dm = pulse_data[:,2]
        self.pulse_sm = pulse_data[:,3]
        self.pulse_spectral_index = pulse_data[:,4]
        self.pulse_data = pulse_data[:,5:]

        self.pulse_sigma = (self.pulse_data - self.noise_mean[np.newaxis,:]) / self.noise_stddev[np.newaxis,:]
        self.pulse_sigmap = (self.pulse_data - self.noise_mean[np.newaxis,:])


    def plot_histogram(self, ialgo, filename, xmax=None):
        import matplotlib.pyplot as plt

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


    def plot_sigma(self, stem, ialgo_list=None, color_list=None, marker_list=None, label_list=None, xmin=None, xmax=None, legloc='lower left'):
        """Writes two files ${stem}_sigma.pdf and ${stem}_sigmap.pdf"""

        import matplotlib.pyplot as plt

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

        todo = [ (stem + '_sigma.pdf', self.pulse_sigma, r'$\sigma$'),
                 (stem + '_sigmap.pdf', self.pulse_sigmap, r"$\sigma'$") ]

        for (filename, d, ylabel) in todo:
            slist = [ ]
            for (i,ialgo) in enumerate(ialgo_list):
                slist.append(plt.scatter(self.pulse_dm, d[:,ialgo], s=5, color=color_list[i], marker=marker_list[i]))

            if xmin is not None:
                plt.xlim(xmin=xmin)
            if xmax is not None:
                plt.xlim(xmax=xmax)

            plt.ylim(ymin=0)
            plt.xlabel('DM')
            plt.ylabel(ylabel)
            plt.legend(slist, label_list, scatterpoints=1, loc=legloc)

            plt.savefig(filename)
            plt.clf()
            print 'wrote', filename
