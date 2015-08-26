#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('hdf5_filename');
parser.add_argument('-o', '--outfile', dest='outfile')

args = parser.parse_args()

if args.outfile:
    import matplotlib
    matplotlib.use('Agg')

import sys
import h5py
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File(args.hdf5_filename, 'r')

data = f['DATA'][:,:]
print data.shape

dmax = float(np.max(data))
dmin = float(np.min(np.where(data >= -1.0e29, data, dmax)))
print 'dmin, dmax =', (dmin,dmax)
assert -1.0e29 < dmin < dmax

mask = np.where(data >= -1.0e29, 1.0, 0.0)
data = (data-dmin) / (dmax-dmin)
data *= mask

colormap = getattr(matplotlib.cm, matplotlib.rcParams['image.cmap'])
rgb = colormap(data)[:,:,:3]
rgb *= mask[:,:,np.newaxis]

plt.imshow(rgb, origin='lower')

if args.outfile:
    plt.savefig(args.outfile)
    print >>sys.stderr, 'wrote', args.outfile
else:
    plt.show()
