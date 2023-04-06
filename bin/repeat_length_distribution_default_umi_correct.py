#!/usr/bin/env python

"""
To plot the UMI distribution pattern.

"""

import sys
import os
from collections import Counter
from statistics import mean
from utils import plot_repeat_dist

csv         = sys.argv[1]
sample_name = sys.argv[2]
outdir      = sys.argv[3]
cutoff      = sys.argv[4]
cutoff      = int(cutoff)
N           = 5
bin_number  = sys.argv[5]

repeat_length_mean = {}
repeat_length_mode = {}
repeat_length_ld   = {} # least distance approach

if os.path.getsize(csv) > 0:
    with open(csv, "r") as f:
        for line in f:
            umi, count, mean, mode, ld = line.split(",")
            count, mean, mode, ld = int(count), float(mean), int(mode), float(ld)

            # support cutoff_0: do not apply any correction
            addon = 1
            if cutoff == 0:
                addon = count
            if count >= cutoff:
                if not mean in repeat_length_mean:
                    repeat_length_mean[mean] = addon
                else:
                    repeat_length_mean[mean] += addon

                if not mode in repeat_length_mode:
                    repeat_length_mode[mode] = addon
                else:
                    repeat_length_mode[mode] += addon

                if not ld in repeat_length_ld:
                    repeat_length_ld[ld] = addon
                else:
                    repeat_length_ld[ld] += addon
else:
    pass
# output to stat:
outfile_mean = os.path.join(outdir, "mean", "stat_mean_" + sample_name + "_cutoff_" + str(cutoff) + ".csv")
with open(outfile_mean, "w") as f:
    for k in sorted(repeat_length_mean.keys()):
        f.write(str(k) + ',' + str(repeat_length_mean[k]) + '\n')
outfile_mode = os.path.join(outdir, "mode", "stat_mode_" + sample_name + "_cutoff_" + str(cutoff) + ".csv")
with open(outfile_mode, "w") as f:
    for k in sorted(repeat_length_mode.keys()):
        f.write(str(k) + ',' + str(repeat_length_mode[k]) + '\n')
outfile_ld = os.path.join(outdir, "ld", "stat_ld_" + sample_name + "_cutoff_" + str(cutoff) + ".csv")
with open(outfile_ld, "w") as f:
    for k in sorted(repeat_length_ld.keys()):
        f.write(str(k) + ',' + str(repeat_length_ld[k]) + '\n')

# output to plot:
outplot_mean = os.path.join(outdir, "mean", "plot_mean_" + sample_name + "_cutoff_" + str(cutoff) + ".png")
if os.path.getsize(outfile_mean) > 0:
    plot_repeat_dist(outfile_mean, outplot_mean, sample_name, N, bin_number)

outplot_mode = os.path.join(outdir, "mode", "plot_mode_" + sample_name + "_cutoff_" + str(cutoff) + ".png")
if os.path.getsize(outfile_mode) > 0:
    plot_repeat_dist(outfile_mode, outplot_mode, sample_name, N, bin_number)

outplot_ld = os.path.join(outdir, "ld", "plot_ld_" + sample_name + "_cutoff_" + str(cutoff) + ".png")
if os.path.getsize(outfile_ld) > 0:
    plot_repeat_dist(outfile_ld, outplot_ld, sample_name, N, bin_number)
