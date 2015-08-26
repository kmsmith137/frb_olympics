#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('algo_filename')
parser.add_argument('output_stem')

args = parser.parse_args()

import frb_olympics
import matplotlib
matplotlib.use('Agg')

data = frb_olympics.compare_run(args.output_stem, algo_filename=args.algo_filename)
data.plot_sigma(args.output_stem)
