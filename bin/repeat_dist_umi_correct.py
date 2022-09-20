#!/usr/bin/env python

"""
Plot UMI distribution.
"""

import sys
import os
from collections import Counter
from statistics import mean
from utils import plot_repeat_dist

tsv         = sys.argv[1]
sample_name = sys.argv[2]
outdir      = sys.argv[3]
cutoff      = sys.argv[4]
cutoff      = int(cutoff)
N           = 5

repeat_length_mean = {}
repeat_length_mode = {}

with open(tsv, "r") as f:
    for line in f:
        umi, count, mean, mode = line.split()
        count, mean, mode = int(float(count)), int(float(mean)), int(mode)
        if count >= cutoff:
            if not mean in repeat_length_mean:
                repeat_length_mean[mean] = 1
            else:
                repeat_length_mean[mean] += 1

            if count >= cutoff:
                if not mode in repeat_length_mode:
                    repeat_length_mode[mode] = 1
                else:
                    repeat_length_mode[mode] += 1
# output to stat:
outfile_mean = os.path.join(outdir, "stat_mean_" + sample_name + "_cutoff_" + str(cutoff) + ".tsv")
with open(outfile_mean, "w") as f:
    for k in sorted(repeat_length_mean.keys()):
        f.write(str(k) + '\t' + str(repeat_length_mean[k]) + '\n')
outfile_mode = os.path.join(outdir, "stat_mode_" + sample_name + "_cutoff_" + str(cutoff) + ".tsv")
with open(outfile_mode, "w") as f:
    for k in sorted(repeat_length_mode.keys()):
        f.write(str(k) + '\t' + str(repeat_length_mode[k]) + '\n')

# output to plot:
outplot_mean = os.path.join(outdir, "plot_mean_" + sample_name + "_cutoff_" + str(cutoff) + ".png")
plot_repeat_dist(outfile_mean, outplot_mean, sample_name, N)

outplot_mode = os.path.join(outdir, "plot_mode_" + sample_name + "_cutoff_" + str(cutoff) + ".png")
plot_repeat_dist(outfile_mode, outplot_mode, sample_name, N)
