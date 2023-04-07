#!/usr/bin/env python

"""
Dev notes:
1.  pd.cut follows: (] range by default

"""

import sys
import os
import numpy as np
import pandas as pd

allele_number  = int(sys.argv[1])
sample_name    = sys.argv[2]
repeat_length_count_table = sys.argv[3]
outfile        = sys.argv[4]
start_allele_1 = sys.argv[5]
end_allele_1   = sys.argv[6]
start_allele_2 = sys.argv[7]
end_allele_2   = sys.argv[8]

# 1. read in sample column into df and add bin:
df = pd.read_csv(repeat_length_count_table)[["repeat_length", sample_name]]
if allele_number == 1:
    bins = [min(df["repeat_length"]) - 1, int(start_allele_1), int(end_allele_1), max(df["repeat_length"]) + 1]
elif allele_number == 2:
    bins = [min(df["repeat_length"]) - 1, int(start_allele_1), int(end_allele_1), int(start_allele_2), int(end_allele_2), max(df["repeat_length"]) + 1]
df["bin"] = pd.cut(df["repeat_length"], bins = bins)

# 2. calculate count, fraction, mean, std:
_total_count = df[sample_name].sum()
df["weighted_length"] = df["repeat_length"] * df[sample_name]

_count = df.groupby("bin")[sample_name].sum()
_fraction = _count / _total_count
_mean = df.groupby("bin")["weighted_length"].sum() / df.groupby("bin")[sample_name].sum()
_std = df.groupby("bin")["repeat_length"].agg(lambda x: np.sqrt(np.average((x - np.average(x, weights = df.loc[x.index, sample_name])) ** 2, weights = df.loc[x.index, sample_name])))

# 3. output to csv:
res = ""
if allele_number == 1:
    header = "below,allele_1,above_allele_1"
elif allele_number == 2:
    header = "below,allele_1,between,allele_2,above"
for i in range(len(_count)):
    # res = res + ",".join([str(_count.iloc[i]), str(_fraction.iloc[i]), str(_mean.iloc[i]), str(_std.iloc[i])])
    res = res + str(_count.iloc[i])
    if not i == len(_count) - 1:
        res = res + ","

with open(outfile, "w") as f:
    f.write(header + "\n")
    f.write(res)
