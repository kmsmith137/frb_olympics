import os
import sys
import numpy as np

from frb_olympics_c import frb_rng, frb_pulse, frb_search_params, \
    simple_direct, simple_tree, sloth, bonsai


algo_list = [ ]
memhack_list = [ ]


def add_algo(algo, memhack=1):
    required_fields = [ 'name', 'debug_buffer_ndm', 'debug_buffer_nt', 'search_gb', 
                        'search_start', 'search_chunk', 'search_end' ]

    missing_fields = [ f for f in required_fields if not hasattr(algo,f) ]

    if len(missing_fields) > 0:
        raise RuntimeError('algorithm object is missing the following required fields: %s' % missing_fields)

    assert isinstance(algo.name, basestring)
    assert len(algo.name) > 0
    assert algo.debug_buffer_ndm > 0
    assert algo.debug_buffer_nt > 0
    assert algo.search_gb >= 0.0
    assert memhack >= 1

    algo_list.append(algo)
    memhack_list.append(memhack)


class compare_run:
    def __init__(self, output_stem, algo_filename=None):
        noise_filename = output_stem + '_noise.txt'
        self.noise_data = np.loadtxt(noise_filename)

        if self.noise_data.ndim == 1:
            self.noise_data = np.reshape(self.noise_data, (-1,1))

        assert self.noise_data.ndim == 2

        self.nnoise = self.noise_data.shape[0]
        self.nalgo = self.noise_data.shape[1]
        self.noise_mean = np.mean(self.noise_data, axis=0)
        self.noise_stddev = np.std(self.noise_data, axis=0)

        pulse_filename = output_stem + '_pulse.txt'
        pulse_data = np.loadtxt(pulse_filename)
        assert pulse_data.shape == (pulse_data.shape[0], self.nalgo+5)

        self.npulse = pulse_data.shape[0]
        self.pulse_arrival_time = pulse_data[:,0]
        self.pulse_intrinsic_width = pulse_data[:,1]
        self.pulse_dm = pulse_data[:,2]
        self.pulse_sm = pulse_data[:,3]
        self.pulse_spectral_index = pulse_data[:,4]
        self.pulse_data = pulse_data[:,5:]

        self.pulse_sigma = (self.pulse_data - self.noise_mean[np.newaxis,:]) / self.noise_stddev[np.newaxis,:]
        self.pulse_sigmap = (self.pulse_data - self.noise_mean[np.newaxis,:])

        if algo_filename is not None:
            self.algo_name_list = [ ]

            for line in open(algo_filename):
                i = line.find('#')
                if i >= 0:
                    line = line[0:i]
    
                line = line.split()
                if len(line) > 0:
                    line = '-'.join(line)
                    line = line.replace('_','-')
                    self.algo_name_list.append(line)
                    
            assert len(self.algo_name_list) == self.nalgo


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
        print >>sys.stderr, 'wrote', filename


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
            assert hasattr(self, 'algo_name_list')
            label_list = [ self.algo_name_list[ialgo] for ialgo in ialgo_list ]

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
            print >>sys.stderr, 'wrote', filename


####################################################################################################
#
# Utility routines


def imp(filename):
    import imp
    module_name = os.path.basename(filename)
    module_name = module_name[:module_name.find('.')]
    return imp.load_source(module_name, filename)


