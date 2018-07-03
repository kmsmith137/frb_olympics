This is a "production-scale" CHIME example.  To run it, you will need a fair amount of memory (say 8-10 GB).
The example assumes that you're using one of the CHIME nodes, say `frb1.physics.mcgill.ca`, which has the
CHIME bonsai configuration files in `/data/bonsai_configs`.

To run the example with 4 Monte Carlo simulations:
```
run-frb-olympics -N -o example3 -n 4 search_params.json bonsai_production_ups_nbeta2_v2.json
```
This won't be enough simulations to get useful statistics, but is a starting point for larger runs.
It will take a while to run!

Some notes:

  - If you look at the parameter file `bonsai_production_ups_nbeta2_v2.json` in this directory,
    you'll see that it is just a pointer to the bonsai config file `/data/bonsai_configs/bonsai_production_ups_nbeta2_v2.hdf5`.
    This is one of four "production-scale" config files that we use in the real search.  The four config files differ in whether
    they use one trial spectral index or two (`nbeta1` vs `nbeta2`), and whether they do or do not use an "upsampled tree" to
    improve sensitivity to FRB's with DM < 820 (`ups` vs `noups`). 

  - All of the "production-scale" bonsai configs search to max DM 13000, corresponding to 5-minute dispersion delay
    across the CHIME band.

    For this reason, in `search_params.json`, I set `nsamples=600000`, corresponding to 10 minutes of data.
    However, I decided to set `dm_max=2000`, to make the sims more representative of real FRB's.
