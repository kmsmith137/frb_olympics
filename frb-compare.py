#!/usr/bin/env python

import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sn', type=float, default=30.0)
parser.add_argument('-o', '--outstem')
parser.add_argument('search_params_txtfile')
parser.add_argument('algo_pyfile')
parser.add_argument('nmc_noise', type=int)
parser.add_argument('nmc_pulse', type=int)

args = parser.parse_args()

if args.outstem is not None:
    output_stem = args.outstem
else:
    assert args.algo_pyfile[-3:] == '.py'
    assert len(args.algo_pyfile) >= 4
    output_stem = args.algo_pyfile[:-3]

assert args.sn > 0.0
assert args.nmc_noise >= 0
assert args.nmc_pulse >= 0


####################################################################################################


import time
import numpy as np
import frb_olympics

frb_olympics.init_mpi_log_files(output_stem)

print 'command line: %s' % (' '.join(sys.argv))
print 'mpi_rank=%d, mpi_size=%d' % (frb_olympics.mpi_rank, frb_olympics.mpi_size)
print 'mpi_rank_within_node=%d, mpi_tasks_per_node=%d' % (frb_olympics.mpi_rank_within_node, frb_olympics.mpi_tasks_per_node)

assert args.nmc_noise % frb_olympics.mpi_size == 0
assert args.nmc_pulse % frb_olympics.mpi_size == 0
(nmc_noise_tot, nmc_pulse_tot) = (args.nmc_noise, args.nmc_pulse)

nmc_noise_loc = nmc_noise_tot // frb_olympics.mpi_size
nmc_pulse_loc = nmc_pulse_tot // frb_olympics.mpi_size


search_params = frb_olympics.frb_search_params(args.search_params_txtfile)
frb_olympics.imp(args.algo_pyfile)
frb_olympics.init_algorithms(search_params)
nalgo = len(frb_olympics.algo_list)

print '--------------------  Search params  --------------------'
search_params.write()

print '--------------------  Memory analysis  --------------------'

timestream_gb = 1.0e-9 * search_params.nchan * search_params.nsamples_per_chunk * 4.0
print 'timestream chunk size = %s GB' % timestream_gb

search_gb = 0.0
for algo in frb_olympics.algo_list:
    print '    %s %s GB' % (algo.name, algo.search_gb)
    search_gb = max(search_gb, algo.search_gb)

print 'total size = %s GB' % (timestream_gb + search_gb)

chunk = np.zeros((search_params.nchan,search_params.nsamples_per_chunk), dtype=np.float32)
score_noise = np.zeros((nmc_noise_loc,nalgo), dtype=np.float)
score_pulse = np.zeros((nmc_pulse_loc,nalgo), dtype=np.float)
pulse_tab = np.zeros((nmc_pulse_loc,5), dtype=np.float)
rng = frb_olympics.frb_rng()


for inoise in xrange(nmc_noise_loc):
    print 'starting noise sim %d/%d' % (inoise, nmc_noise_loc)
    rsave = frb_olympics.frb_rng(rng)

    for (ialgo, algo) in frb_olympics.enumerate_algorithms_with_memhack():
        rng = frb_olympics.frb_rng(rsave)
        t0 = time.time()

        for ichunk in xrange(search_params.nchunks):
            print '    %s: starting chunk %s/%s' % (algo.name, ichunk, search_params.nchunks)
            search_params.simulate_noise(rng, chunk)
            algo.search_chunk(chunk, ichunk)

        print '    %s: search_result=%s [time=%s]' % (algo.name, algo.search_result, time.time()-t0)
        score_noise[inoise,ialgo] = algo.search_result


for ipulse in xrange(nmc_pulse_loc):
    print 'starting pulse sim %d/%d' % (ipulse, nmc_pulse_loc)

    p = search_params.make_random_pulse(rng, 1.0)
    p.fluence = args.sn / search_params.get_signal_to_noise_of_pulse(p)

    pulse_tab[ipulse,0] = p.arrival_time
    pulse_tab[ipulse,1] = p.intrinsic_width
    pulse_tab[ipulse,2] = p.dispersion_measure
    pulse_tab[ipulse,3] = p.scattering_measure
    pulse_tab[ipulse,4] = p.spectral_index

    for (ialgo, algo) in frb_olympics.enumerate_algorithms_with_memhack():
        t0 = time.time()

        for ichunk in xrange(search_params.nchunks):
            print '    %s: starting chunk %s/%s' % (algo.name, ichunk, search_params.nchunks)
            chunk[:,:] = 0.0
            search_params.add_pulse(p, chunk, ichunk)
            algo.search_chunk(chunk, ichunk)

        print '    %s: search_result=%s [time=%s]' % (algo.name, algo.search_result, time.time()-t0)
        score_pulse[ipulse,ialgo] = algo.search_result


# lists of numpy arrays
score_noise = frb_olympics.mpi_gather(score_noise)
score_pulse = frb_olympics.mpi_gather(score_pulse)
pulse_tab = frb_olympics.mpi_gather(pulse_tab)


if frb_olympics.mpi_rank == 0:
    noise_data = np.concatenate(score_noise)
    noise_filename = args.outstem + '_noise.txt'
    np.savetxt(noise_filename, noise_data)
    print 'wrote %s' % noise_filename

    pulse_data = np.concatenate((np.concatenate(pulse_tab), np.concatenate(score_pulse)), axis=1)
    pulse_filename = args.outstem + '_pulse.txt'

    f = open(pulse_filename, 'w')
    print >>f, '# col 0: arrival time'
    print >>f, '# col 1: intrinsic width'
    print >>f, '# col 2: dispersion measure'
    print >>f, '# col 3: scattering measure'
    print >>f, '# col 4: spectral index'

    np.savetxt(f, pulse_data)
    print 'wrote %s' % pulse_filename

    import matplotlib
    matplotlib.use('Agg')

    c = frb_olympics.comparison_outputs(noise_data, pulse_data)
    c.plot_sigma(args.outstem)
