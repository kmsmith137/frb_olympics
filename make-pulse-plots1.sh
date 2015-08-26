#!/bin/bash
#
# For visual testing of pulse generation code
#
# Reminder: write-pulse command-line syntax is
#   write-pulse <arrival_time> <intrinsic_width> <DM> <SM> <freq_lo> <freq_hi> <dt_sample> <nsamples_per_chunk>
#
# All pulses have fluence 1 and specral index 0.


#
# Case 1: realistically dispersed pulse, plotted at two frequencies and a variety of sampling rates
#
# Note: dispersion delay at DM=600 and nu=700 MHz is 5.080 sec
#
mkdir -p pulse1

./write-pulse 10.0 0.0 600.0 0.0 700 701 1.0e-3 4 > pulse1/freq1_dt1.txt
./write-pulse 10.0 0.0 600.0 0.0 700 701 2.038271e-3 2 > pulse1/freq1_dt2.txt
./write-pulse 10.0 0.0 600.0 0.0 700 701 3.017312e-3 3 > pulse1/freq1_dt3.txt
./write-pulse 10.0 0.0 600.0 0.0 700 701 4.029831e-3 1 > pulse1/freq1_dt4.txt
./write-pulse 10.0 0.0 600.0 0.0 700 701 5.058189e-3 2 > pulse1/freq1_dt5.txt

./write-pulse 10.0 0.0 600.0 0.0 703 704 1.0e-3 1 > pulse1/freq2_dt1.txt
./write-pulse 10.0 0.0 600.0 0.0 703 704 2.0e-3 3 > pulse1/freq2_dt2.txt
./write-pulse 10.0 0.0 600.0 0.0 703 704 3.0e-3 2 > pulse1/freq2_dt3.txt
./write-pulse 10.0 0.0 600.0 0.0 703 704 4.0e-3 2 > pulse1/freq2_dt4.txt
./write-pulse 10.0 0.0 600.0 0.0 703 704 5.0e-3 3 > pulse1/freq2_dt5.txt


#
# Case 2: scattered pulses with no dispersion
#

mkdir -p pulse2
./write-pulse 8.0 0.0 0.0 10.0 700 701 1.0e-3 40 > pulse2/sm10.txt
./write-pulse 8.0 0.0 0.0 13.0 700 701 1.151e-3 31 > pulse2/sm13.txt
./write-pulse 8.0 0.0 0.0 16.0 700 701 1.429e-3 55 > pulse2/sm16.txt
./write-pulse 8.0 0.0 0.0 19.0 700 701 1.681e-3 67 > pulse2/sm19.txt

./write-pulse 8.3 0.05 0.0 0.0 700 701 0.7321e-3 33 > pulse2/wsm0.txt
./write-pulse 8.3 0.05 0.0 5.0 700 701 1.0483e-3 49 > pulse2/wsm5.txt
./write-pulse 8.3 0.05 0.0 10.0 700 701 1.238e-3 37 > pulse2/wsm10.txt
./write-pulse 8.3 0.05 0.0 15.0 700 701 1.581e-3 64 > pulse2/wsm15.txt
