#!/usr/bin/env python

import frb_olympics

r = frb_olympics.olympics('example0_search_params.json', simulate_noise=False)
r.add_bonsai('example0_bonsai_config_n2048.hdf5', name='2K-trees')
r.add_bonsai('example0_bonsai_config_n4096.hdf5', name='4K-trees')

# run() writes the output files 'example0.json' and 'example0_snr_vs_dm.pdf'
r.run('example0.json', 16)
