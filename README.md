### frb_olympics

### Dependencies

The core `frb_olympics` module has the following dependencies:

   - numpy
   - matplotlib
   - [kmsmith137/simpulse](https://github.com/kmsmith137/simpulse).  Note that simpulse requires cython and FFTW3.

You'll also want one or more dedispersers.  These are optional dependencies: you can run `frb_olympics`
on whatever subset of them you have installed.

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
