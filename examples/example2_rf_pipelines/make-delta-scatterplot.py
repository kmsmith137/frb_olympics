#!/usr/bin/env python
#
# This script postprocesses the frb_olympics output file (example2.json) and produces a
# "non-standard" scatterplot with SNR on the x-axis,and Delta on the y-axis, where Delta is defined by:
#
#    Delta = DM / (pulse_width**0.5 + sample_width^2)**0.5


import frb_olympics

e = frb_olympics.ensemble.from_json('example2.json')

e.make_snr_plot('Delta.pdf',
                lambda d: d['dm'] / (1.0e6*d['intrinsic_width']**2 + 1)**0.5,
                r'$\Delta$',
                xmin = 0,
                xmax = 80,
                legloc = 'lower right')
