#!/usr/bin/env python

"""
Calculate the fraction of reads given four cutoffs.

sample_name, below_count, below_frac, below_mean, below_sd, between_count, between_frac, between_sd, above_count, above_frac, above_mean, above_sd, below_2_count, below_2_frac, below_2_mean, below_2_sd, between_2_count, between_2_frac, between_2_sd, above_2_count, above_2_frac, above_2_mean, above_2_sd

(below_count) cutoff_1 (bewteen_count) cutoff_2 (above_count | below_2_count) cutoff_3 (between_2_count) cutoff_4 (above_2_count)

"""

import sys
import os
import numpy as np

sample_name = sys.argv[1]
csv = sys.argv[2]
outdir = sys.argv[3]
umi_cutoff = sys.argv[4]
cutoff_below, cutoff_above, cutoff_below_2, cutoff_above_2 = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])

count_below, count_between, count_above, count_between_2, count_above_2, count_all = 0, 0, 0, 0, 0, 0
below_list, between_list, above_list, bewteen_list_2, above_list_2 = [], [], [], [], []

if os.path.getsize(csv) > 0:
    with open(csv, "r") as f:
        for line in f:
            length, c = line.strip().split(",")
            if not length in ["plus", "problem"]:
                length, c = int(float(length)), int(float(c))
                count_all = count_all + c
                if length < cutoff_below:
                    count_below += c
                    below_list = below_list + [length] * c
                if length >= cutoff_below and length < cutoff_above:
                    count_between += c
                    between_list = between_list + [length] * c
                if length >= cutoff_above and length < cutoff_below_2:
                    count_above += c
                    count_above = above_list + [length] * c
                if length >= cutoff_below_2 and length < cutoff_above_2:
                    count_between_2 += c
                    between_list_2 = between_list_2 + [length] * c
                if length >= cutoff_above_2:
                    count_above_2 += c
                    above_list_2 = above_list_2 + [length] * c
            elif length == "plus":
                count_all += int(c)
                count_above += int(c)
else:
    below_count, below_frac, below_mean, below_std = "nan", "nan", "nan", "nan"
    above_count, above_frac, above_mean, above_std = "nan", "nan", "nan", "nan"
    between_count, between_frac, between_mean, between_std = "nan", "nan", "nan", "nan"

outfile = os.path.join(outdir, sample_name + "_frac_" + str(cutoff_below) + "_" + str(cutoff_above) + "_cutoff_" + umi_cutoff + ".csv")

with open(outfile, "w") as f:
    if os.path.getsize(csv) > 0:
        below_count = str(count_below)
        below_frac  = str(count_below/count_all)
        below_mean  = str(np.mean(below_list))
        below_std   = str(np.std(below_list))
        between_count = str(count_between)
        between_frac  = str(count_between/count_all)
        between_mean  = str(np.mean(between_list))
        between_std   = str(np.std(between_list))
        above_count = str(count_above)
        above_frac  = str(count_above/count_all)
        above_mean  = str(np.mean(above_list))
        above_std   = str(np.std(above_list))

    res_list = [sample_name, below_count, below_frac, below_mean, below_std, between_count, between_frac, between_mean, between_std, above_count, above_frac, above_mean, above_std]
    res = ",".join(res_list) + "\n"
    f.write(res)
