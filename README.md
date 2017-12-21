### frb_olympics

A simulation framework for studying optimality of FRB detection algorithms.

The `frb_olympics` code is pure python, but it calls code written in low-level languages (e.g. bonsai).

### Dependencies

The core `frb_olympics` module has the following dependencies:

   - numpy
   - matplotlib
   - [kmsmith137/simpulse](https://github.com/kmsmith137/simpulse).  Note that simpulse requires cython and FFTW3.

You'll also want one or more dedispersers.  These are optional dependencies: you can run `frb_olympics`
using whatever subset of them you have installed.

   - bonsai (https://github.com/CHIMEFRB/bonsai).  Note that bonsai has its own dependencies:
      - cython
      - [kmsmith137/simd_helpers](https://github.com/kmsmith137/simd_helpers):
        header-only library for writing x86 assembly language kernels.
      - [kmsmith137/rf_kernels](https://github.com/kmsmith137/rf_kernels):
        fast C++/assembly kernels for RFI removal and related tasks.

   - bz_fdmt (https://github.com/kmsmith137/bz_fdmt).

### Installation:

This is a standard `setuptools` install:
```
python setup.py install --prefix=$HOME
```
To verify that it worked:

  - Switch to a different directory, and do `python -c 'import frb_olympics'` to check that python
    interpreter is finding the frb_olympics module.  If not, you need to add the setuptools install
    directory (see setuptools output, should be something like `$HOME/lib/python2.7/site-packages`)
    to your `$PYTHONPATH`.

  - Do `run-frb-olympics` to check that the shell is finding the run-frb-olympics script.  If not,
    you need to add the setuptools install directory (see setuptools output, should be something like
    `$HOME/bin`) to your `$PATH`.

### Documentation

For a quick start, see examples in the `examples/` directory.

Detailed documentation is in the python docstrings.

Some high-level points to be aware of:

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
    **milliseconds** (not seconds).  This is the only place where we use milliseconds instead
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
