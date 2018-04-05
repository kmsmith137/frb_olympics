#!/usr/bin/env python
#
# Creates 'rfp_bonsai.json' and 'rfp_detrender.json'.
#
# These are rf_pipelines json files (not frb_olympics json files)
# which represent serialized rf_pipelines.pipeline_objects.

import rf_pipelines

d = rf_pipelines.bonsai_dedisperser(config_filename = 'bonsai_ntree4096_nups2.hdf5',
                                    fill_rfi_mask = False,
                                    use_analytic_normalization = True,
                                    track_global_max = True)

rf_pipelines.json_write('rfp_bonsai.json', d, clobber=True)

d1 = rf_pipelines.polynomial_detrender(nt_chunk=1024, axis='time', polydeg=4)
d2 = rf_pipelines.spline_detrender(nt_chunk=1024, axis='freq', nbins=6)
p = rf_pipelines.pipeline([d1,d2])

rf_pipelines.json_write('rfp_detrender.json', p, clobber=True)
