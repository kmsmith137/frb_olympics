#!/usr/bin/env python

import frb_olympics_x

r = frb_olympics_x.olympics('example0_search_params.txt', simulate_noise=False)
r.add_bonsai('example0_bonsai_config.hdf5', name='example0')
r.run('example0.json', 128)
