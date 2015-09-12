#!/usr/bin/env python
#
# This is a throwaway script from when I was visually debugging the pulse simulation code.
# I may take it out of git, or replace it with a more rigorous test.

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import frb_olympics


def plot_pulse(arrival_time, freq_lo, freq_hi, dt_sample, nsamples_per_chunk, intrinsic_width=0.0, dm=0.0, sm=0.0):
    p = frb_olympics.frb_pulse(1.0, arrival_time, intrinsic_width, dm, sm, 0.0)   # fluence=1, spectral_index=0

    (t0,t1) = p.get_endpoints(freq_lo, freq_hi)
    i0 = int(t0/dt_sample) - 2
    i1 = int(t1/dt_sample) + 3
    
    nchunks = int(i1/nsamples_per_chunk) + 2
    ts = np.zeros(nchunks * nsamples_per_chunk, dtype=np.float32)
    
    for ichunk in xrange(nchunks):
        i = (ichunk) * nsamples_per_chunk
        j = (ichunk+1) * nsamples_per_chunk
        p.add_to_timestream(freq_lo, freq_hi, ts[i:j], dt_sample, ichunk)

    x = np.zeros(2*(i1-i0))
    y = np.zeros(2*(i1-i0))
    x[0::2] = dt_sample * np.arange(i0,i1,dtype=np.float)
    x[1::2] = dt_sample * np.arange(i0+1,i1+1,dtype=np.float)
    y[0::2] = ts[i0:i1]
    y[1::2] = ts[i0:i1]

    plt.plot(x,y)


#
# Case 1: realistically dispersed pulse, plotted at two frequencies and a variety of sampling rates
#
# Note: dispersion delay at DM=600 and nu=700 MHz is 5.080 sec
#

plot_pulse(10., 700., 701., 1.0e-3, 4, dm=600.)
plot_pulse(10., 700., 701., 2.038271e-3, 2, dm=600.)
plot_pulse(10., 700., 701., 3.017312e-3, 3, dm=600.)
plot_pulse(10., 700., 701., 4.029831e-3, 1, dm=600.)
plot_pulse(10., 700., 701., 5.058189e-3, 2, dm=600.)

plot_pulse(10., 703., 704., 1.0e-3, 1, dm=600.)
plot_pulse(10., 703., 704., 2.0e-3, 3, dm=600.)
plot_pulse(10., 703., 704., 3.0e-3, 2, dm=600.)
plot_pulse(10., 703., 704., 4.0e-3, 2, dm=600.)
plot_pulse(10., 703., 704., 5.0e-3, 3, dm=600.)

plt.savefig('pulse1.pdf')
plt.close()
print 'wrote pulse1.pdf'

#
# Case 2: scattered pulses with no dispersion
#

plot_pulse(8.0, 700., 701., 1.0e-3, 40, sm=10.0)
plot_pulse(8.0, 700., 701., 1.151e-3, 31, sm=13.0)
plot_pulse(8.0, 700., 701., 1.429e-3, 55, sm=16.0)
plot_pulse(8.0, 700., 701., 1.681e-3, 67, sm=19.0)

# plot_pulse(8.3, 701., 701., 0.7321e-3, 33, sm=0.01)
plot_pulse(8.3, 701., 701., 1.0483e-3, 49, sm=5.0)
plot_pulse(8.3, 701., 701., 1.238e-3, 37, sm=10.0)
plot_pulse(8.3, 701., 701., 1.581e-3, 64, sm=15.0)

plt.savefig('pulse2.pdf')
plt.close()
print 'wrote pulse2.pdf'
