### frb_olympics

There used to be a lot of code here, but all functionality was moved to other repositories.
All that's left is a 500-line python file!

### Dependencies

Strictly speaking, the only dependencies are the following libraries:

   - simpulse (https://github.com/kmsmith137/simpulse)
   - rf_pipelines (https://github.com/kmsmith137/rf_pipelines)
   - bonsai (https://github.com/CHIMEFRB/bonsai)

However, the complete list of dependencies for the above libraries reads as follows:

   - C++ compiler which supports C++11

   - A very recent cython.  I know that cython 0.24 works, and cython 0.20 is too old.
     A hint for upgrading cython: `[sudo] pip upgrade [--user] Cython`


   - FFTW3 (http://fftw.org)

   - libhdf5 (https://www.hdfgroup.org/HDF5/release/obtain5.html)
     Note that this is a link to HDF5 v1.8.  I imagine v1.10 also works but haven't tested it yet.

   - jsoncpp (https://github.com/open-source-parsers/jsoncpp)

     On osx, you can probably install with: `brew install jsoncpp`
 
     In linux, you can probably install with `yum install jsoncpp-devel`

     Building jsoncpp from scratch is a pain, but the following procedure worked for me:
     ```
     git clone https://github.com/open-source-parsers/jsoncpp
     mkdir -p build/debug
     cd build/debug
     cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_C_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=debug -G "Unix Makefiles" ../..
     make install
     ```

   - Optional but recommended: The 'PIL' python imaging library (you can test whether you have 
     it with 'import PIL' from python).  If you need to install it, I recommend the 'Pillow' 
     variant (pip install Pillow)

### Installation:

To install frb_olympics, just create a Makefile.local defining the variable PYDIR where rf_pipelines.py should be installed.
(See examples in site/)  Then do 'make install'.

From there, try running examples/example0_bonsai.  (The other examples are intended to run on a large cluster and are more complicated.)

