#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outfile', dest='outfile')
parser.add_argument('-i', '--ialgo', type=int, default=0)
parser.add_argument('search_params_txtfile')
parser.add_argument('algo_pyfile')

args = parser.parse_args()

if args.outfile:
    import matplotlib
    matplotlib.use('Agg')

import sys
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import frb_olympics

search_params = frb_olympics.frb_search_params(args.search_params_txtfile)
frb_olympics.imp(args.algo_pyfile)
frb_olympics.init_algorithms(search_params)

assert 0 <= args.ialgo < len(frb_olympics.algo_list)
algo = frb_olympics.algo_list[args.ialgo]

print 'Debug buffer shape = (%d,%d)' % (algo.debug_buffer_ndm, algo.debug_buffer_nt)

gb = algo.search_gb
gb += 1.0e-9 * search_params.nchan * search_params.nsamples_per_chunk * 4.0   # timestream
gb += 1.0e-9 * algo.debug_buffer_ndm * algo.debug_buffer_nt * 4.0   # timestream
print 'Estimated memory usage: %s GB' % gb

sm = search_params.sm_max
if sm > 0.0:
    print 'Using SM =', sm

(dm0, dm1) = (search_params.dm_min, search_params.dm_max)
dm_list = [ (0.97*dm0 + 0.03*dm1), (0.50*dm0 + 0.50*dm1), (0.03*dm0 + 0.97*dm1) ]

print 'DM range in search_params =', ((dm0,dm1))
print 'Putting pulses at DM list =', dm_list

pulse_list = [ ]
for dm in dm_list: 
    (t0, t1) = search_params.get_allowed_arrival_times(0.0, dm, sm)
    pulse_list.append(frb_olympics.frb_pulse(1.0, 0.97*t0 + 0.03*t1, 1.0e-5, dm, sm, 0.0))
    pulse_list.append(frb_olympics.frb_pulse(1.0, 0.50*t0 + 0.50*t1, 1.0e-5, dm, sm, 0.0))
    pulse_list.append(frb_olympics.frb_pulse(1.0, 0.03*t0 + 0.97*t1, 1.0e-5, dm, sm, 0.0))

debug_buffer = np.zeros((algo.debug_buffer_ndm, algo.debug_buffer_nt), dtype=np.float32)
debug_buffer[:] = -1.0e30

algo.search_start(0)

for ichunk in xrange(search_params.nchunks):
    timestream = np.zeros((search_params.nchan, search_params.nsamples_per_chunk), dtype=np.float32)
    for p in pulse_list:
        search_params.add_pulse(p, timestream, ichunk)

    print '%s: processing chunk %s/%s' % (algo.name, ichunk, search_params.nchunks)
    algo.search_chunk(timestream, ichunk, debug_buffer)
    del timestream

algo.search_end()

dmax = float(np.max(debug_buffer))
dmin = float(np.min(np.where(debug_buffer >= -1.0e29, debug_buffer, dmax)))
print 'dmin, dmax =', (dmin,dmax)
assert -1.0e29 < dmin < dmax

mask = np.where(debug_buffer >= -1.0e29, 1.0, 0.0)
debug_buffer = (debug_buffer-dmin) / (dmax-dmin)
debug_buffer *= mask

colormap = getattr(matplotlib.cm, matplotlib.rcParams['image.cmap'])
rgb = colormap(debug_buffer)[:,:,:3]
rgb *= mask[:,:,np.newaxis]

plt.imshow(rgb, origin='lower')

if args.outfile:
    plt.savefig(args.outfile)
    print >>sys.stderr, 'wrote', args.outfile
else:
    plt.show()
