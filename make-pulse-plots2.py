#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def plot_file(filename):
    assert filename[-4:]=='.txt'
    #h = {'1':'k-', '2':'b-', '3':'r-', '4':'m-', '5':'g-'}
    #fmt = h[filename[-5]]
    a = np.loadtxt(filename)
    n = len(a)
    x = np.zeros(2*n)
    y = np.zeros(2*n)
    x[0::2] = a[:,0]
    x[1::2] = a[:,1]
    y[0::2] = a[:,2]
    y[1::2] = a[:,2]
    plt.plot(x,y)

plot_file('pulse1/freq1_dt1.txt')
plot_file('pulse1/freq1_dt2.txt')
plot_file('pulse1/freq1_dt3.txt')
plot_file('pulse1/freq1_dt4.txt')
plot_file('pulse1/freq1_dt5.txt')

plot_file('pulse1/freq2_dt1.txt')
plot_file('pulse1/freq2_dt2.txt')
plot_file('pulse1/freq2_dt3.txt')
plot_file('pulse1/freq2_dt4.txt')
plot_file('pulse1/freq2_dt5.txt')

plt.show()
plt.close()

plot_file('pulse2/sm10.txt')
plot_file('pulse2/sm13.txt')
plot_file('pulse2/sm16.txt')
plot_file('pulse2/sm19.txt')

plot_file('pulse2/wsm0.txt')
plot_file('pulse2/wsm5.txt')
plot_file('pulse2/wsm10.txt')
plot_file('pulse2/wsm15.txt')

plt.show()
plt.close()
