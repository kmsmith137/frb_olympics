### frb_olympics/example1

To run the example:
```
run-frb-olympics -N -o example1 -n 64 search_params.json bonsai_ntree2048_nups1.json bonsai_ntree4096_nups2.json bonsai_ntree8192_nups4.json bz_fdmt.py
```
The outputs are `example1.json` and `example1_snr_vs_dm.pdf`.

Note that the DM's and arrival times reported by bonsai are usually coarse-grained for speed, and 
therefore contain errors due to coarse-graining.  To disable the coarse-graining, you can set
`dm_coarse_graining` and `time_coarse_graining` to 1 in the bonsai config file.
