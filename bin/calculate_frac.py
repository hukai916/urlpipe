#!/usr/bin/env python

"""
Calculate the fraction of reads given two cutoffs.

sample_name, below_count, below_frac, below_mean, below_sd, above_count, above_frac, above_mean, above_sd

"""

import sys
import os
import numpy as np
import pandas as pd

sample_name = sys.argv[1]
csv = sys.argv[2]
outdir = sys.argv[3]
cutoff_below, cutoff_above = int(sys.argv[4]), int(sys.argv[5])

if os.path.getsize(csv) > 0:
    # csv = "30"
    df = pd.read_csv(csv, header = None, names = ["length", "frequency"])
    # skip "plus" and "problem" reads because they make it tricky to calculate the mean length
    df["length"] = pd.to_numeric(df["length"], errors = 'coerce')
    df.dropna(subset=["length"], inplace = True)
    bins = [0, cutoff_below, cutoff_above, np.inf]
    df["bin"] = pd.cut(df["length"], bins=bins)

    # calculate count, fraction, mean, std:
    _total_count = df["frequency"].sum()
    df["weighted_length"] = df["length"] * df["frequency"]

    _count = df.groupby("bin")["frequency"].sum()
    _fraction = _count / _total_count
    _mean = df.groupby("bin")["weighted_length"].sum() / df.groupby("bin")["frequency"].sum()
    _std = df.groupby("bin")["length"].agg(lambda x: np.sqrt(np.average((x - np.average(x, weights = df.loc[x.index, "frequency"])) ** 2, weights = df.loc[x.index, "frequency"])))

    res = ""
    for i in range(len(_count)):
        res = res + ",".join([str(_count.iloc[i]), str(_fraction.iloc[i]), str(_mean.iloc[i]), str(_std.iloc[i])])
        if not i == len(_count) - 1:
            res = res + ","
else:
    res = ",".join(["nan"] * 12)

outfile = os.path.join(outdir, sample_name + "_frac_" + str(cutoff_below) + "_" + str(cutoff_above) + "_cutoff_" + str(cutoff_below) + "_" + str(cutoff_below) + ".csv")
with open(outfile, "w") as f:
    f.write(res)
    # print(res)
