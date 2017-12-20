import os
import sys
import copy
import json
import importlib

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import simpulse


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
                 dm_min=0.0, sm_min=0.0, sm_max=0.0, spectral_index_min=0.0, spectral_index_max=0.0,
                 intrinsic_width_min=0.0, intrinsic_width_max=0.0, snr_min=20.0, snr_max=20.0, filename=None):

        self.nfreq = nfreq
        self.freq_lo_MHz = freq_lo_MHz
        self.freq_hi_MHz = freq_hi_MHz
        self.nsamples = nsamples
        self.dt_sample = dt_sample
        self.dm_min = dm_min
        self.dm_max = dm_max
        self.sm_min = sm_min
        self.sm_max = sm_max
        self.spectral_index_min = spectral_index_min
        self.spectral_index_max = spectral_index_max
        self.intrinsic_width_min = intrinsic_width_min
        self.intrinsic_width_max = intrinsic_width_max
        self.snr_min = snr_min
        self.snr_max = snr_max

        # Argument checking follows.
        
        if filename is None:
            filename = 'frb_olympics.search_params constructor'

        assert self.dm_min >= 0.0, filename + ": failed assert: dm_min >= 0.0"
        assert self.sm_min >= 0.0, filename + ": failed assert: sm_min >= 0.0"
        assert self.intrinsic_width_min >= 0.0, filename + ": failed assert: intrinsic_width_min >= 0.0"
        assert self.nsamples > 0, filename + ": failed assert: nsamples > 0"
        assert self.nfreq > 0, filename + ": failed assert: nfreq > 0"

        # The choice of ranges here is intended to guard against accidentally using the wrong units
        # (e.g. GHz instead of MHz, millseconds instead of seconds)
        assert self.freq_lo_MHz >= 100.0, filename + ": failed assert: freq_lo_MHz >= 100.0"
        assert self.dt_sample >= 2.0e-6, filename + ": failed assert: dt_sample >= 2.0e-6"
        assert self.dt_sample <= 0.01, filename + ": failed assert: dt_sample <= 0.01"

        assert self.dm_min <= self.dm_max, filename + ": failed assert: dm_min <= dm_max"
        assert self.sm_min <= self.sm_max, filename + ": failed assert: sm_min <= sm_max"
        assert self.intrinsic_width_min <= self.intrinsic_width_max, filename + ": failed assert: intrinsic_width_min <= intrinsic_width_max"
        assert self.freq_lo_MHz < self.freq_hi_MHz, filename + ": failed assert: freq_lo_MHz < freq_hi_MHz"
        assert self.snr_min <= self.snr_max, ": failed assert: snr_min <= snr_max"
        assert self.snr_min > 0, ": failed assert: snr_min > 0"

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
                 'spectral_index_min': self.spectral_index_min,
                 'spectral_index_max': self.spectral_index_max,
                 'intrinsic_width_min': self.intrinsic_width_min,
                 'intrinsic_width_max': self.intrinsic_width_max,
                 'snr_min': self.snr_min,
                 'snr_max': self.snr_max }

    
    @staticmethod
    def from_json(j, filename=None):
        if filename is None:
            filename = 'frb_olympics.search_params.from_json()'
        if not isinstance(j, dict):
            raise RuntimeError('%s: expected dict, got %s' % (filename, j.__class__.__name__))

        required_keys = set(['nfreq', 'freq_lo_MHz', 'freq_hi_MHz', 'nsamples', 'dt_sample', 'dm_max'])
        optional_keys = set(['dm_min', 'sm_min', 'sm_max', 'spectral_index_min', 'spectral_index_max', 'intrinsic_width_min', 'intrinsic_width_max', 'snr_min', 'snr_max' ])
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
    def __init__(self, name):
        self.name = name


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

    
    def jsonize(self):
        raise RuntimeError("frb_olympics.dedisperser_base.jsonize() was not overloaded by subclass %s" % self.__class__.__name__)


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


class comparison:
    def __init__(self, sp, dedisperser_list):
        assert isinstance(sp, search_params)
        assert len(dedisperser_list) > 0

        self.search_params = sp
        self.dedisperser_list = dedisperser_list
        self.dedisperser_json = [ ]
        self.sim_json = [ ]

        for d in dedisperser_list:
            assert isinstance(d, dedisperser_base)
            assert hasattr(d, 'name')

            j = d.jsonize()

            if not isinstance(j, dict):
                raise RuntimeError("expected %s.jsonize() to return dict (returned %s)" % (p.__class__.__name__, j.__class__.__name__))

            j['module_name'] = d.__module__
            j['class_name'] = d.__class__.__name__
            j['name'] = d.name

            self.dedisperser_json.append(j)
            d.init_search_params(sp)


    def run(self, nmc, noisy=True):
        nfreq = self.search_params.nfreq
        nsamples = self.search_params.nsamples
        freq_lo_MHz = self.search_params.freq_lo_MHz
        freq_hi_MHz = self.search_params.freq_hi_MHz
        dm_max = self.search_params.dm_max

        nmc_in = len(self.sim_json)
        intensity = np.zeros((nfreq, nsamples), dtype=np.float32)

        timestream_length = nsamples * self.search_params.dt_sample
        max_dispersion_delay = dispersion_delay(dm_max, freq_lo_MHz) - dispersion_delay(dm_max, freq_hi_MHz)
        noise_seeds = [ ]

        for (id,d) in enumerate(self.dedisperser_list):
            d.allocate()

            for imc in xrange(nmc):
                if noisy:
                    print 'frb_olympics: dedisperser %d/%d (%s), simulation %d' % (id+1, len(self.dedisperser_list), d.name, nmc_in+imc+1)

                if id == 0:
                    # Simulate random FRB params.
                    true_params = {
                        'dm': np.random.uniform(self.search_params.dm_min, self.search_params.dm_max),
                        'sm': np.random.uniform(self.search_params.sm_min, self.search_params.sm_max),
                        'spectral_index': np.random.uniform(self.search_params.spectral_index_min, self.search_params.spectral_index_max),
                        'intrinsic_width': np.random.uniform(self.search_params.intrinsic_width_min, self.search_params.intrinsic_width_max),
                        'snr': np.random.uniform(self.search_params.snr_min, self.search_params.snr_max),
                        'tmid': np.random.uniform(0.51 * max_dispersion_delay, timestream_length - 0.51 * max_dispersion_delay)
                    }

                    this_sim = {
                        'true_params': true_params,
                        'recovered_params': [ ]    # this list will be populated by the dedispersers
                    }

                    # Save FRB params and RNG state
                    self.sim_json.append(this_sim)
                    noise_seeds.append(np.random.get_state())

                else:
                    # Use same random FRB params and RNG state as previous simulation.
                    true_params = self.sim_json[nmc_in + imc]['true_params']
                    np.random.set_state(noise_seeds[imc])

                # Simulate Gaussian random noise.  There is no float32 gaussian random number generator in numpy, 
                # so we simulate in float64 and down-convert.  We do this in slices to save memory!

                for i in xrange(nsamples):
                    intensity[:,i] = np.random.standard_normal(size=nfreq)

                # Add simulated FRB.  Note that we pay the computational cost of simulating the pulse
                # "from scratch" for every dedisperser.  This sometimes (if nmc is large) saves memory,
                # and the cost of simulating the pulse is small.

                # The simpulse library uses the undispersed arrival time 'tu' of the pulse,
                # whereas frb_olympics uses the central arrival time 'tmid', so we need to translate.

                dt_i = dispersion_delay(true_params['dm'], freq_hi_MHz)
                dt_f = dispersion_delay(true_params['dm'], freq_lo_MHz)
                true_tu = true_params['tmid'] - dt_i - (dt_f-dt_i)/2.0

                # Note that 'nt' is the number of samples used internally by simpulse to represent
                # the pulse.  (FIXME: tune this parameter.) 

                p = simpulse.single_pulse(nt = 1024,
                                          nfreq = nfreq,
                                          freq_lo_MHz = freq_lo_MHz,
                                          freq_hi_MHz = freq_hi_MHz,
                                          dm = true_params['dm'],
                                          sm = true_params['sm'],
                                          intrinsic_width = true_params['intrinsic_width'],
                                          fluence = 1.0,
                                          spectral_index = true_params['spectral_index'],
                                          undispersed_arrival_time = true_tu)

                # We simulate the pulse with an nominal normalization (fluence=1.0), and
                # rescale to the target SNR.

                nominal_snr = p.get_signal_to_noise(self.search_params.dt_sample)
                rescaling_factor = true_params['snr'] / nominal_snr
                p.fluence *= rescaling_factor

                # Add FRB here!
                p.add_to_timestream(intensity, 0.0, timestream_length, freq_hi_to_lo=True)

                # Run dedisperser here!
                dedisperser_output = d.dedisperse(intensity)

                # The return value of d.dedisperser() should be a dictionary containing 'snr', 'dm',
                # and precisely one of { 'tmid', 'tini', or 'tfin' }.

                if not isinstance(dedisperser_output, dict):
                    raise RuntimeError('expected return value of %s.dedisperse() to be a dict, got %s' %  (d.__class__.__name__, dedisperser_output.__class__.__name__))

                required_keys = set(['snr','dm'])
                optional_keys = set(['tini','tmid','tfin'])
                missing_keys = required_keys.difference(dedisperser_output.keys())
                unrecognized_keys = set(dedisperser_output.keys()).difference(required_keys).difference(optional_keys)

                if len(missing_keys) > 0:
                    raise RuntimeError('return value of %s.dedisperse() does not contain key(s) %s' %  (d.__class__.__name__, sorted(missing_keys)))
                if len(unrecognized_keys) > 0:
                    raise RuntimeError('return value of %s.dedisperse() contains unrecognized key(s) %s' %  (d.__class__.__name__, sorted(unrecognized_keys)))
                if len(dedisperser_output.keys()) != 3:
                    raise RuntimeError('return value of %s.dedisperse() must contain precisely one of %s' % (d.__class__.__name__, sorted(optional_keys)))

                # If tini or tfin were specified, translate to a value of tmid.

                dm = dedisperser_output['dm']
                dt = dispersion_delay(dm, freq_lo_MHz) - dispersion_delay(dm, freq_hi_MHz)

                if not dedisperser_output.has_key('tmid') and dedisperser_output.has_key('tini'):
                    dedisperser_output['tmid'] = dedisperser_output['tini'] + dt/2.0
                if not dedisperser_output.has_key('tmid') and dedisperser_output.has_key('tfin'):
                    dedisperser_output['tmid'] = dedisperser_output['tfin'] - dt/2.0

                assert dedisperser_output.has_key('tmid')
                
                # Done with this (simulation, dedisperser) pair.
                # Record the results in the 'sim_json' data structure.

                self.sim_json[nmc_in+imc]['recovered_params'].append(dedisperser_output)

            d.deallocate()


    def jsonize(self):
        return {
            'search_params': self.search_params.jsonize(),
            'dedisperser_list': self.dedisperser_json,
            'sims': self.sim_json
        }


    @staticmethod
    def from_json(j):
        sp = search_params.from_json(j['search_params'])
        dlist = [ dedisperser_base.from_json(j) for j in j['dedisperser_list'] ]
        return comparison(sp, dlist)


    def make_snr_plot(self, plot_filename, xaxis_param, xaxis_label, legend_labels = None):
        if len(self.sim_json) <= 1:
            print '%s: no plot written, not enough sims' % plot_filename
            return

        xmin = getattr(self.search_params, xaxis_param + '_min')
        xmax = getattr(self.search_params, xaxis_param + '_max')

        if xmin == xmax:
            print "%s: no plot written, the parameter '%s' was not varied in this run" % (plot_filename, xaxis_param)
            return

        if (legend_labels == None) or (legend_labels == [ ]):
            legend_labels = [ ]

            for d in self.dedisperser_list:
                # Initialize legend_labels from dedisperser names.
                # Dedisperser names often contain underscores, which confuse matplotlib's tex rendering,
                # so we replace '_' by r'\_'.  This is not a systematic approach to making the names
                # "TeX safe", but it's usually enough in practice!

                n = copy.copy(d.name)
                n = n.replace('_', r'\_')
                legend_labels.append(n)

        xvec = np.array([ s['true_params'][xaxis_param] for s in self.sim_json ])
        yarr = np.array([ [ (r['snr']/s['true_params']['snr']) for r in s['recovered_params'] ] for s in self.sim_json ])

        assert xvec.shape == (len(self.sim_json),)
        assert yarr.shape == (len(self.sim_json), len(self.dedisperser_list))
        assert len(legend_labels) == len(self.dedisperser_list)

        plt.xlim(xmin, xmax)
        plt.ylim(0.0, max(np.max(yarr),1.1))
        plt.axhline(y=1, ls=':', color='k')
        plt.xlabel(xaxis_label)
        plt.ylabel('Optimality')

        slist = [ ]
        color_list = [ 'b', 'r', 'm', 'g', 'k', 'y', 'c']

        for ids in xrange(len(self.dedisperser_list)):
            c = color_list[ids % len(color_list)]
            s = plt.scatter(xvec, yarr[:,ids], s=5, color=c, marker='o')
            slist.append(s)

            # Dedisperser names often contain underscores, which confuse matplotlib's tex rendering,
            # so we replace '_' by r'\_'.  This is not a systematic approach to making the names
            # "TeX safe", but it's usually enough in practice!

            n = self.dedisperser_list[ids].name
            n = n.replace('_', r'\_')
            legend_labels.append(n)


        plt.legend(slist, legend_labels, scatterpoints=1, loc='lower left')
        plt.savefig(plot_filename)
        plt.clf()

        print 'wrote', plot_filename
        

    def make_snr_plots(self, plot_filename_stem, legend_labels = None):
        todo = [ ('dm', 'DM'),
                 ('sm', 'SM'),
                 ('spectral_index', 'Spectral index'),
                 ('intrinsic_width', 'Intrinsic width') ]

        for (xaxis_param, xaxis_label) in todo:
            plot_filename = '%s_snr_vs_%s.pdf' % (plot_filename_stem, xaxis_param)
            self.make_snr_plot(plot_filename, xaxis_param, xaxis_label, legend_labels)
