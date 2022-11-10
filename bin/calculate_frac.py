#!/usr/bin/env python

"""
Calculate the fraction of reads given two cutoffs.

sample_name, below_count, below_frac, below_mean, below_sd, above_count, above_frac, above_mean, above_sd

"""

import sys
import os
import numpy as np

sample_name = sys.argv[1]
csv = sys.argv[2]
outdir = sys.argv[3]
umi_cutoff = sys.argv[4]
cutoff = sys.argv[5] # below/above cutoff

count_below, count_above, count_all = 0, 0, 0
cutoff_below, cutoff_above = cutoff.strip().split()
cutoff_below, cutoff_above = int(cutoff_below), int(cutoff_above)
below_list, above_list = [], []

with open(csv, "r") as f:
    for line in f:
        length, c = line.strip().split(",")
        if not length in ["plus", "problem"]:
            length, c = int(float(length)), int(float(c))
            count_all = count_all + c
            if length < cutoff_below:
                count_below += c
                below_list = below_list + [length] * c
            if length > cutoff_above:
                count_above += c
                above_list = above_list + [length] * c
        elif length == "plus":
            count_all += int(c)
            count_above += int(c)

outfile = os.path.join(outdir, sample_name + "_frac_" + str(cutoff_below) + "_" + str(cutoff_above) + "_cutoff_" + umi_cutoff + ".csv")
with open(outfile, "w") as f:
    below_count = str(count_below)
    below_frac  = str(count_below/count_all)
    below_mean  = str(np.mean(below_list))
    below_std   = str(np.std(below_list))
    above_count = str(count_above)
    above_frac  = str(count_above/count_all)
    above_mean  = str(np.mean(above_list))
    above_std   = str(np.std(above_list))

    res_list = [sample_name, below_count, below_frac, below_mean, below_std, above_count, above_frac, above_mean, above_std]
    res = ",".join(res_list) + "\n"
    f.write(res)
