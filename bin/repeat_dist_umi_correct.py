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

repeat_length_mean = {}
repeat_length_mode = {}

with open(tsv, "r") as f:
    for line in f:
        umi, count, mean, mode = line.split()
        if count >= int(cutoff):
            if not int(mean) in repeat_length_mean:
                repeat_length_mean[int(mean)] = 1
            else:
                repeat_length_mean[int(mean)] += 1

            if count >= int(cutoff):
                if not int(mode) in repeat_length_mode:
                    repeat_length_mode[int(mode)] = 1
                else:
                    repeat_length_mode[int(mode)] += 1
# output to stat:
with open(os.path.join(outdir, "stat_mean_" + smaple_name + "_cutoff_" + str(cutoff) + ".tsv"), "w") as f:
    for k in sorted(repeat_length_mean.keys()):
        f.write(str(k) + '\t' + str(repeat_length_mean[k]) + '\n')
with open(os.path.join(outdir, "stat_mode_" + smaple_name + "_cutoff_" + str(cutoff) + ".tsv"), "w") as f:
    for k in sorted(repeat_length_mode.keys()):
        f.write(str(k) + '\t' + str(repeat_length_mode[k]) + '\n')


with open(outdir + "/test.txt", "w") as f:
    f.write("Simply a test!")
    f.write(print_test())
