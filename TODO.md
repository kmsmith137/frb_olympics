Here is a "wish list" of dedispersers we'd like to compare:

  - Heimdall (in progress)

  - sigproc

  - presto `single_pulse_search`

  - FREDDA (not public)

  - AstroAccelerate

General:

  - Currently, the entire intensity array is held in memory, and dedispersed with a single call
    dedisperser.dedisperse().  For large timestreams, it would be better to simulate the array
    incrementally, and define an API for incremental dedispersion.

  - Currently, we implement Gaussian profiles (if the `intrinsic_width` parameter is nonzero).
    Is this the best choice?  Would it make sense to have a boolean flag to switch to boxcar profiles?

  - Multiprocessing/MPI runs

  - Command-line tools for combining/appending 

  - More plotting and postprocessing tools (for example, based on success fraction rather than
    reported SNR.)

Bonsai:

  - Currently, the analytic transfer matrix is computed from scratch whenever a bonsai_dedisperser
    is constructed, which is annoying!  Should cache it in an HDF5 file (this is mostly implemented
    in bonsai already.)

  - Cleanup: less verbose output from bonsai_dedisperser.jsonize()
