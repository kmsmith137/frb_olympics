"""
frb_olympics: Simulation framework for studying optimality of FRB detection algorithms.

For a quick start, try the examples in the `examples/` directory.

Detailed documentation is in the python docstrings.  Some high-level points to be aware of:

  - An frb_olympics Monte Carlo ensemble (`class ensemble`) consists of:

      - A set of parameter ranges (represented by `class search_params`) for the 
        following parameters: DM, pulse width, spectral index, scattering timescale,
        signal-to-noise ratio.  In each Monte Carlo realization, each parameter is
        randomly sampled from its allowed range.

      - A list of dedispersers that will run on the simulations.  These are subclasses
        of `class dedisperser_base`, and are "thin" wrappers around code in other
        repositories.  

        Currently, we implement `class bonsai_dedipserser`, which is a wrapper around 
        bonsai (https://github.com/CHIMEFRB/bonsai), and `class bz_fdmt_dedisperser`,
        which is a wrapper around FDMT (https://github.com/kmsmith137/bz_fdmt).
        More coming soon!

      - In each Monte Carlo realization, a simulation is generated consisting of
        Gaussian random noise, plus a random FRB.  Each dedisperser analyzes the
        simulation, makes one guess for the location of the FRB, and returns three
        numbers: the estimated DM, arrival time, and SNR.

      - The driver script `run-frb-olympics` will run an ensemble of Monte Carlo
        simulations, and write a JSON file containing the search_params, the
        dedisperser list, the true FRB parameters in each simulation, and the
        recovered parameters from each dedisperser.

      - This JSON file can be postprocessed to produce various plots.

  - Scattering is implemented as an exponential profile whose characteristic timescale
    depends on frequency as f^(-4.4).  In contrast, the "intrinsic width" of an FRB is
    implemented as a frequency-independent Gaussian.

    In each frequency channel, the pulse shape is the convolution of these two profiles,
    plus a boxcar profile which represents dispersion delay within the channel.
    
    We should decide whether the Gaussian intrinsic profile is the best choice.  For
    example, we could 

  - We define the scattering measure (SM) to be the scattering timescale at 1 GHz, in
    MILLISECONDS (not seconds).  This is the only place where we use milliseconds instead
    of seconds!

  - There are four possible definitions of the arrival time of an FRB:

      - "initial" arrival time: arrival time at the highest frequency in the band (i.e. least delayed)
      - "final" arrival time: arrival time at the lowest frequency in the band (i.e. most delayed)
      - "middle" arrival time: average of initial and final times (warning: not the arrival time at the central frequency!)
      - "undispersed" arrival time: arrival time in the limit of high frequency.

    In the core frb_olympics code, we generally use t_middle, but the individual dedisperser
    classes can return either t_initial, t_middle, or t_final, and frb_olympics will translate 
    to a value of t_middle.

  - The 'run-frb-olympics' script has a -N flag which deserves special discussion.  If
    specified, then the simulations will contain an FRB with no noise.  (By default, if -N
    is not specified, then the simulations will contain an FRB + noise.)
    
    This option only produces reasonable results if all of the dedispersers use precomputed
    variances to normalize their signal-to-noise.  This is the case for both of the dedispersers
    currently implemented (`bonsai_dedisperser` and `bz_fdmt_dedisperser`), so using the -N
    flag makes the SNR plots look a little nicer, by removing noise scatter.

    However, for many dedispersers (such as Heimdall) the variances are estimated directly
    from the output of the dedispersion transform, rather than being precomputed.  In this
    case, using the -N flag will result in spuriously large SNR values, and results will not
    make sense!
"""

import os
import sys
import copy
import json
import importlib

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import simpulse


def dispersion_delay(dm, freq_MHz):
    """Returns dispersion delay in seconds, for given DM and frequency."""

    return 4.148806e3 * dm / (freq_MHz * freq_MHz);


class json_read_helper:
    """
    This helper class is used by functions which take a json object argument 'j'
    for example search_params.from_json(j).  It allows 'j' to be either a python
    dictionary or a filename.  It also does some error checking.

    Members:
    
       self.json: python dictionary representing json data

       self.filename: name of file that was originally read to obtain 'self.json'.
          This can be None, if no filename is available (e.g. if json was constructed
          directly from python, rather than reading from a file).

       self.diagnostic_name: a string which can be printed during error reporting.
          This is never None.
    """

    def __init__(self, j, filename, caller_name, expected_keys = [ ]):
        """
        Constructor arguments:

          j: either a python dictionary representing json data, or a string
             which will be interpreted as a json filename.  (Be careful not
             to pass the string serialization of the json data, or it will
             be interpreted as a filename!)

          filename: if 'j' is a dictionary was originally read from a file,
             this should be the name of the file.  Otherwise, it can be None.

          caller_name: name of function which called json_read_helper().
             This cannot be none.
        
          expected_keys: if this is a nonempty list of strings, then some
             error checking will be performed, by checking the json dictionary
             for the presence of the given keys.
        """

        assert caller_name is not None

        if isinstance(j, basestring):
            # If 'j' is a string, interpret as filename
            filename = j
            f = open(filename)
    
            try:
                j = json.load(f)
            except:
                raise RuntimeError("%s: couldn't parse json file" % filename)

        self.json = j
        self.filename = filename   # can be None
        self.diagnostic_name = filename if (filename is not None) else caller_name

        if not isinstance(j, dict):
            raise RuntimeError('%s: expected dict, got %s' % (self.diagnostic_name, j.__class__.__name__))

        for k in expected_keys:
            if not j.has_key(k):
                raise RuntimeError("%s: key '%s' not found" % (self.diagnostic_name, k))


####################################################################################################


class search_params:
    """
    search_params: represents a set of parameters for an frb_olympics run.

       self.nfreq: number of frequency channels (int)
       self.freq_lo_MHz: lowest frequency in band (float, MHz)
       self.freq_hi_MHz: lowest frequency in band (float, MHz)
       self.nsamples: number of time samples in each Monte Carlo simulation (int)
       self.dt_sample: length of a time sample (float, seconds)
       self.dm_min: minimum dispersion measure (DM) used in the frb_olympics run (float, pc cm^(-3))
       self.dm_max: maximum dispersion measure (DM) used in the frb_olympics run (float, pc cm^(-3))
       self.sm_min: minimum scattering measure (SM) used in the frb_olympics run (float, milliseconds)
       self.sm_max: maximum scattering measure (SM) used in the frb_olympics run (float, milliseconds)
       self.spectral_index_min: minimum spectral index used in the frb_olympics run (float, dimensionless)
       self.spectral_index_max: maximum spectral index used in the frb_olympics run (float, dimensionless)
       self.intrinsic_width_min: minimum intrinsic width used in the frb_olympics run (float, seconds)
       self.intrinsic_width_max: maximum intrinsic width used in the frb_olympics run (float, seconds)
       self.snr_min: minimum signal-to-noise ratio used in the frb_olympics run (float, dimensionless)
       self.snr_max: maximum signal-to-noise ratio used in the frb_olympics run (float, dimensionless)

    When an FRB is simulated, its parameters (DM, SM, SNR, intrinsic width, spectral index)
    will be randomly drawn from the ranges above.

    Note that we define the scattering measure (SM) to be the scattering timescale at 1 GHz, 
    in MILLISECONDS (not seconds).  This is the only place where we use milliseconds instead
    of seconds!

    The JSON file format for the seach_params is just a "flat" dictionary, see `examples/*/search_params.json`
    for examples.
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
        """Returns a python dictionary which is a valid JSON object."""

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
        """The "inverse" of jsonize(): it takes a python dictionary which is valid JSON, and returns a search_params instance."""

        r = json_read_helper(j, filename, 'frb_olympics.search_params.from_json()')

        required_keys = set(['nfreq', 'freq_lo_MHz', 'freq_hi_MHz', 'nsamples', 'dt_sample', 'dm_max'])
        optional_keys = set(['dm_min', 'sm_min', 'sm_max', 'spectral_index_min', 'spectral_index_max', 'intrinsic_width_min', 'intrinsic_width_max', 'snr_min', 'snr_max' ])

        actual_keys = set([ str(x) for x in r.json.keys() ])
        missing_keys = required_keys.difference(actual_keys)
        unrecognized_keys = set(actual_keys).difference(required_keys).difference(optional_keys)

        if len(missing_keys) > 0:
            raise RuntimeError('%s: missing key(s): %s' % (r.diagnostic_name, sorted(missing_keys)))
        if len(unrecognized_keys) > 0:
            raise RuntimeError('%s: unrecognized key(s): %s' % (r.diagnostic_name, sorted(unrecognized_keys)))

        kwds = copy.copy(r.json)
        kwds['filename'] = filename

        return search_params(**kwds)


####################################################################################################


class dedisperser_base:
    """
    dedipserser_base: this abstract base class represents a dedisperser.

    Subclasses of dedisperser_base are "thin" wrappers around code in other repositories,
    for example `class bonsai_dedisperser` wraps bonsai (https://github.com/CHIMEFRB/bonsai).

    Subclasses must implement the following functions:

       self.init_search_params(sparams)    required
       self.allocate()                     optional
       self.dedisperse(intensity)          required
       self.deallocate()                   optional
       self.jsonize()                      required
       self.from_json()                    required, staticmethod

    These functions might be called in a sequence which looks something like this:

       d = dedisperser_subclass.from_json()
       d.init_search_params(sparams)    # only called once
       
       # First "block" of simulations
       d.allocate()
       d.dedisperse(intensity_mc1)
       d.dedisperse(intensity_mc2)
       d.deallocate()

       # Second "block" of simulations
       d.allocate()
       d.dedisperse(intensity_mc3)
       d.dedisperse(intensity_mc4)
       d.deallocate()
    """

    def __init__(self, tex_label, precomputed_variance):
        """
        The subclass constructor should call this base class constructor.

        The 'precomputed_variance' flag should be True if the dedisperser uses precomputed variance
        to normalize its signal-to-noise estimates (e.g. bonsai, FDMT).  It should be False if the
        dedisperser estimates its variance directly from the output arrays (e.g. heimdall).
        
        This detail matters because 'run-frb-olympics -N' only makes sense if all dedispersers
        use precomputed variances.  (For more discussion, see 'run-frb-olympics -h'.)
        """

        assert tex_label is not None

        self.tex_label = tex_label
        self.precomputed_variance = precomputed_variance


    def init_search_params(self, sparams):
        """Must be overridden by subclass.  The 'sparams' argument is an instance of 'class search_params'."""
        raise RuntimeError("frb_olympics.dedisperser_base.init_search_params() was not overridden by subclass %s" % self.__class__.__name__)


    def allocate(self):
        """Overriding this in the subclass is optional, but should probably be done if the dedipserser uses a lot of memory or other resources."""
        pass


    def dedisperse(self, intensity):
        """
        Must be overridden by subclass.  The 'intensity' argument is a float32 array of shape (nfreq, nsamples).

        Note that 'nfreq', 'nsamples', ... are members of the search_params object, which is passed to init_search_params()
        before dedipserse() gets called.

        The return value of dedisperse() is a dictionary with 3 members: 'snr', 'dm', and one of { 'tmid', 'tini', 'tfin' }.
        This allows the dedisperser to use one of three possible definitions of the arrival time of an FRB:
           - "initial" arrival time: arrival time at the highest frequency in the band (i.e. least delayed)
           - "final" arrival time: arrival time at the lowest frequency in the band (i.e. most delayed)
           - "middle" arrival time: average of initial and final times (warning: not the arrival time at the central frequency!)

        (Note that the caller uses 'tmid', and will translate whatever arrival time is returned to a value of 'tmid',
        but this detail shouldn't affect the implementation of dedisperse() in the subclass.)
        """
        raise RuntimeError("frb_olympics.dedisperser_base.dedisperse() was not overridden by subclass %s" % self.__class__.__name__)


    def deallocate(self):
        """Overriding this in the subclass is optional, but should probably be done if the dedipserser uses a lot of memory or other resources."""
        pass

    
    def jsonize(self):
        """
        Must be overridden by subclass.  The return value should be a python dictionary which is valid JSON.

        (Note that the caller will add additional members 'module_name', 'class_name', 'tex_label' to the dictionary,
        but this detail shouldn't affect the implementation of jsonize() in the subclass.)
        """
        raise RuntimeError("frb_olympics.dedisperser_base.jsonize() was not overridden by subclass %s" % self.__class__.__name__)


    @staticmethod
    def from_json(j, filename=None):
        """
        Must be defined by subclass.  This is the "inverse" of jsonize(): it takes a python dictionary which
        is valid JSON, and returns an instance of the subclass.

        One possible point of confusion: there will be two versions of from_json(), the base class staticmethod
        dedisperser_base.from_json(), and the subclass staticmethod dedisperser_subclass.from_json().  To construct
        a dedisperser from json data, the caller first calls dedisperser_base.from_json().  This function determines
        the correct subclass (using the 'module_name' and 'class_name' json members) and then calls the appropriate
        dedisperser_subclass.from_json().

        Assuming this usage, the implementation of from_json() in the subclass can assume that j is a dictionary
        which defines the key 'tex_label' (in addition to 'module_name' and 'class_name', but these probably won't
        be useful to the subclass).
        """

        expected_keys = [ 'module_name', 'class_name' ]
        r = json_read_helper(j, filename, 'frb_olympics.dedisperser_base.from_json()', expected_keys)

        j = copy.copy(r.json)
        module_name = j['module_name']
        class_name = j['class_name']

        # If the tex_label is not specified in the json data, construct one from either
        # the filename (first choice) or class name (second choice).

        if not j.has_key('tex_label'):
            if r.filename is not None:
                t = os.path.basename(r.filename)
                t = t.split('.')[0]
            else:
                t = class_name

            # Replace underscores by r'\_', to avoid crashing latex!
            t = t.replace('_', r'\_')
            j['tex_label'] = t

        # Using the JSON fields 'module_name' and 'class_name', import the appropriate
        # dedisperser_base subclass.  There is a lot of sanity checking here!
            
        try:
            m = importlib.import_module(module_name)
        except ImportError:
            raise ImportError("%s: couldn't import module %s" % (r.diagnostic_name, module_name))

        c = getattr(m, class_name, None)

        if c is None:
            raise RuntimeError("%s: couldn't find class '%s' in module '%s" % (r.diagnostic_name, class_name, module_name))
        if not issubclass(c, dedisperser_base):
            raise RuntimeError("%s: expected class %s.%s to be a subclass of frb_olympics.dedisperser_base" % (r.diagnostic_name, module_name, class_name))

        # This sanity check is important to prevent an infinite loop, if the subclass does not
        # define the from_json() staticmethod.

        if c.from_json == dedisperser_base.from_json:
            raise RuntimeError("%s: expected class %s.%s to override dedisperser_base.from_json()" % (r.diagnostic_name, module_name, class_name))

        # Now delegate to the subclass from_json() staticmethod.
        
        return c.from_json(j, r.filename)


####################################################################################################


class ensemble:
    """
    An frb_olympics Monte Carlo ensemble (`class ensemble`) consists of:

      - A set of parameter ranges (represented by `class search_params`) for the 
        following parameters: DM, pulse width, spectral index, scattering timescale,
        signal-to-noise ratio.  In each Monte Carlo realization, each parameter is
        randomly sampled from its allowed range.

      - A list of dedispersers that will run on the simulations.  These are subclasses
        of `class dedisperser_base`, and are "thin" wrappers around code in other
        repositories.  

        Currently, we implement `class bonsai_dedipserser`, which is a wrapper around 
        bonsai (https://github.com/CHIMEFRB/bonsai), and `class bz_fdmt_dedisperser`,
        which is a wrapper around FDMT (https://github.com/kmsmith137/bz_fdmt).
        More coming soon!

      - In each Monte Carlo realization, a simulation is generated consisting of
        Gaussian random noise, plus a random FRB.  Each dedisperser analyzes the
        simulation, makes one guess for the location of the FRB, and returns three
        numbers: the estimated DM, arrival time, and SNR.
    """

    def __init__(self, sparams, dedisperser_list, sim_json=None, add_noise=True):
        """
        Constructor arguments:

          - sparams: instance of 'class search_params'.
          - dedisperser_list: list of instances of subclasses of 'class dedisperser_base'.
          - add_noise: if True then simulations will be (FRB+noise), if False then sims will be FRB-only.
          - sim_json (optional): json data for previously-run sims, if any.
        """

        assert isinstance(sparams, search_params)
        assert len(dedisperser_list) > 0

        # FIXME should error-check 'sim_json'

        self.search_params = sparams
        self.dedisperser_list = dedisperser_list
        self.dedisperser_json = [ ]
        self.sim_json = sim_json if (sim_json is not None) else [ ]
        self.add_noise = add_noise

        for d in dedisperser_list:
            assert isinstance(d, dedisperser_base)
            
            for k in [ 'tex_label', 'precomputed_variance' ]:
                if not hasattr(d, k):
                    raise RuntimeError("%s: no '%s' member found, you probably forgot to call the base class constructor dedisperser_base.__init__()" % (d.__class__.__name__, k))

            if (not add_noise) and (not d.precomputed_variance):
                raise RuntimeError("%s: this dedisperser does not use precomputed variances, but is being used in a noiseless run (for more info, see 'run-frb-olympics -h')" % d.__class__.__name___)

            j = d.jsonize()

            if not isinstance(j, dict):
                raise RuntimeError("expected %s.jsonize() to return dict (returned %s)" % (p.__class__.__name__, j.__class__.__name__))

            j['module_name'] = d.__module__
            j['class_name'] = d.__class__.__name__
            j['tex_label'] = d.tex_label

            self.dedisperser_json.append(j)

            # Note call dedisperser_base.init_search_params() here.
            d.init_search_params(sparams)


    def run(self, nmc, noisy=True):
        """
        Runs 'nmc' Monte Carlo simulations, and appends the results to the 'class ensemble' internal state.

        In principle the interface is general enough to make either the simulation or the dedisperser the
        outer loop.  For example:
        
            # Run 1000 Monte Carlo simulations, with outer loop over simulation and inner loop over dedisperser.
            ensemble.run(1000)

            # Run 1000 Monte Carlo simulations, with outer loop over dedipserser and inner loop over simulation.
            for i in xrange(1000):
                ensemble.run(1)

        However, I ended up deciding to use the latter form (outer loop over dedipserser and inner loop over simulation)
        exclusively, so the interface could be simplified by replacing run(nmc) by run_one().
        """

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

        # Outer loop over dedispersers.
        for (id,d) in enumerate(self.dedisperser_list):
            d.allocate()

            # Inner loop over simulations.
            for imc in xrange(nmc):
                if noisy:
                    print 'frb_olympics: dedisperser %d/%d (%s), simulation %d' % (id+1, len(self.dedisperser_list), d.tex_label, nmc_in+imc+1)

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

                if self.add_noise:
                    for i in xrange(nsamples):
                        intensity[:,i] = np.random.standard_normal(size=nfreq)
                else:
                    intensity[:,:] = 0.

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
                # Note that we order frequency channels from lowest to highest (the intuitive ordering,
                # but some dedispersers assume the opposite and will need to reverse it, for example
                # bonsai and heimdall).
                p.add_to_timestream(intensity, 0.0, timestream_length)

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

            # End of inner loop over simulations; back to outer loop over dedispersers.
            d.deallocate()


    def make_snr_plot(self, plot_filename, xaxis_param, xaxis_label, legend_labels = None):
        """
        Makes one SNR plot, with 

          - "Optimality" on the y-axis, i.e. (SNR)_recovered / (SNR)_true
          - On the x-axis, either DM, SM, SNR_true, spectral index, or intrinsic width.
          
        Arguments:
        
          - plot_filename: should end in '.pdf'
          - xaxis_param: either 'dm', 'sm', 'snr', 'spectral_index', 'intrinsic_width'
          - xaxis_label: a TeX axis label consistent with xaxis_param, for example 'DM'.
          - legend_labels: list of tex labels, one for each dedisperser.
              (optional: if unspecified, then default tex_labels will be used.)
        """

        if len(self.sim_json) <= 1:
            print '%s: no plot written, not enough sims' % plot_filename
            return

        if legend_labels is None:
            legend_labels = [ d.tex_label for d in self.dedisperser_list ]
            
        assert len(legend_labels) == len(self.dedisperser_list)

        xmin = getattr(self.search_params, xaxis_param + '_min')
        xmax = getattr(self.search_params, xaxis_param + '_max')

        if xmin == xmax:
            print "%s: no plot written, the parameter '%s' was not varied in this run" % (plot_filename, xaxis_param)
            return
                               
        xvec = np.array([ s['true_params'][xaxis_param] for s in self.sim_json ])
        yarr = np.array([ [ (r['snr']/s['true_params']['snr']) for r in s['recovered_params'] ] for s in self.sim_json ])

        assert xvec.shape == (len(self.sim_json),)
        assert yarr.shape == (len(self.sim_json), len(self.dedisperser_list))

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

        plt.legend(slist, legend_labels, scatterpoints=1, loc='lower left')
        plt.savefig(plot_filename)
        plt.clf()

        print 'wrote', plot_filename
        

    def make_snr_plots(self, plot_filename_stem, legend_labels = None):
        """
        Makes all possible SNR plots, with

          - "Optimality" on the y-axis, i.e. (SNR)_recovered / (SNR)_true
          - On the x-axis, either DM, SM, SNR_true, spectral index, or intrinsic width.
          
        Arguments:
        
          - plot_filename_stem: this will be used to generate filenames such as
              ${plot_filename_stem}_snr_vs_${quantity}.pdf
            where quantity = dm, sm, ...

          - legend_labels: list of tex labels, one for each dedisperser.
              (optional: if unspecified, then default tex_labels will be used.)

        Note: if a parameter is not actually varied in the simulations, then the
        corresponding plot will not be written.
        """

        todo = [ ('dm', 'DM'),
                 ('sm', 'SM'),
                 ('spectral_index', 'Spectral index'),
                 ('intrinsic_width', 'Intrinsic width') ]

        for (xaxis_param, xaxis_label) in todo:
            plot_filename = '%s_snr_vs_%s.pdf' % (plot_filename_stem, xaxis_param)
            self.make_snr_plot(plot_filename, xaxis_param, xaxis_label, legend_labels)


    def jsonize(self):
        """Returns a python dictionary which is a valid JSON object."""

        return {
            'add_noise': self.add_noise,
            'search_params': self.search_params.jsonize(),
            'dedisperser_list': self.dedisperser_json,
            'sims': self.sim_json
        }


    @staticmethod
    def from_json(j, filename=None):
        """The "inverse" of jsonize(): it takes a python dictionary which is valid JSON, and returns a search_params instance."""

        expected_keys = [ 'add_noise', 'search_params', 'dedisperser_list', 'sims' ]
        r = json_read_helper(j, filename, 'frb_olympics.ensemble.from_json()', expected_keys)

        sparams = search_params.from_json(r.json['search_params'], r.filename)
        sim_json = r.json['sims']

        dlist = [ ]
        for dj in r.json['dedisperser_list']:
            d = dedisperser_base.from_json(dj, r.filename)
            dlist.append(d)

        return ensemble(sparams, dlist, sim_json)
