This directory contains an example of `rf_pipelines_dedisperser`.

We compare recovered SNR from (1) bonsai alone, (2) detrender + bonsai.
This was part of the `18-03-19-pulsar-snr` study, sent on the CHIME email list on
2008 March 19.

Note that we use nfreq=1024 and dm_max=200, rather than CHIME parameters (nfreq=16384, dm_max=13000),
so that the simulations will run faster.

This is conceptually straightforward, but keeping all the different config files 
straight requires some explanation, since the bonsai, rf_pipelines, and frb_olympics
packages each define different file formats!  To help keep this straight, we use a 
convention where filenames always begin with 'bonsai', 'rfp', or 'foly'.

```
# The input file 'bonsai_ntree4096_nups2.txt' specifies parameters for the bonsai dedisperser.
# The 'bonsai-mkweight' utility is part of bonsai, and creates the bonsai HDF5 file.

bonsai-mkweight bonsai_ntree4096_nups2.txt bonsai_ntree4096_nups2.hdf5


# The 'make-rfp-json-files.py' script creates:
#   - rfp_bonsai.json       (dedisperser object)
#   - rfp_detrender.json    (detrender object)
#
# These are rf_pipelines json files (not frb_olympics json files) which represent serialized
# objects of class rf_pipelines.pipeline_object.  For more details, see the 'make-rfp-json-files.py'
# script itself.

./make-rfp-json-files.py


# Now we create frb_olympics json files, which represent serialized objects of class
# frb_olympics.dedisperser_base.   (In this case, they will be objects of class
# frb_olympics.rf_pipelines_dedisperser, a subclass of frb_olympics.dedisperser_base
# which "wraps" an rf_pipelines pipeline_object).
#
# We use the 'rfp-to-frb-olympics' script, which is part of frb_olympics.  It takes one
# or more rf_pipelines json files, constructs an object of class frb_olympics.rf_pipelines_dedisperser,
# and serializes to a json file.  For more information, do `rfp-to-frb-olympics -h`.

rfp-to-frb-olympics -f -t no_detrender -o foly_no_detrender.json rfp_bonsai.json 
rfp-to-frb-olympics -f -t with_detrender -o foly_with_detrender.json rfp_detrender.json rfp_bonsai.json 


# Now we run the FRB olympics.  Note that 'foly_search_params.json' is an frb_olympics configuration file
# which specifies the parameter ranges in the simulation.
#
# The 'run-frb-olympics' script is part of frb_olympics, and produces the output file example2.json,
# along with a few standard scatterplots.

run-frb-olympics -o example2 -N -n 256 foly_search_params.json foly_no_detrender.json foly_with_detrender.json


# One more minor detail.  In the CHIME study, I wanted to make a non-standard scatterplot with SNR on the x-axis,
# and Delta on the y-axis, where Delta is defined by:
#
#   Delta = DM / (pulse_width**0.5 + sample_width^2)**0.5

./make-delta-scatterplot.py
```
