### frb_olympics/example1

This is the same as example 0, except that we now include FDMT.
In particular, we use 1024 frequency channels from 400-800 MHz, and max DM 200.
Each Monte Carlo simulation contains 8 seconds of data with 1 ms sampling.

To run the example:
```
run-frb-olympics -N -o example1 -n 64 search_params.json bonsai_ntree2048_nups1.json bonsai_ntree4096_nups2.json bonsai_ntree8192_nups4.json bz_fdmt.json
```
The outputs are `example1.json` and `example1_snr_vs_dm.pdf`.

Some notes:

  - The DM's and arrival times reported by bonsai are usually coarse-grained for speed, and 
    therefore contain errors due to coarse-graining.  To disable the coarse-graining, you can set
   `dm_coarse_graining` and `time_coarse_graining` to 1 in the bonsai config file.

  - The tex labels are specified in the bonsai json files (along with the usual bonsai config_params,
    which are explained in the bonsai documentation)

  - We call 'run-frb-olympics' with the -N flag, which deserves special discussion.  If
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
