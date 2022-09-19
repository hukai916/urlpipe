#!/usr/bin/env python

"""
To visualize the repeat length distribution by UMI group.
"""

import sys
import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import pandas as pd

tsv = sys.argv[1]
sample_name = sys.argv[2]
cutoff      = sys.argv[3]
outdir_stat = sys.argv[4]
output_plot = sys.argv[5]
group_num   = sys.argv[6]

# read tsv into dict
d_umi = {}
for line in open(tsv):
    umi, length = line.split()
    umi = umi.split("_")[1]
    if not umi in d_umi:
        d_umi[umi] = [length]
    else:
        d_umi[umi].append(length)

d_cutoff = {}
for x in d_umi:
    if len(d_umi[x]) == int(cutoff):
        if not len(d_cutoff) > int(group_num):
            d_cutoff[x] = d_umi[x]
        else:
            break

# save to stat
output_stat = os.path.join(outdir_stat, sample_name + ".tsv")
os.makedirs(os.path.dirname(output_stat), exist_ok=True)
with open(output_stat, "w") as f:
    for k in d_cutoff:
        for c in d_cutoff[k]:
            f.write(k + "\t" + str(c) + "\n")

# save to plot
df = pd.read_csv(output_stat, sep = "\t", header = None)
for group_id, group_df in df.groupby(df.iloc[:,0]):
    print(group_df.iloc[:,1].tolist())
    plt.hist(group_df.iloc[:,1], bins = max(int(group_df.iloc[:,1].tolist())) - min(int(group_df.iloc[:,1].tolist())))

    plt.xlabel("repeat length")
    plt.ylabel("count")
    plt.title("Repeat Length Distribution within UMI Group: " + group_df.iloc[0,0])

    _output_plot = os.path.join(output_plot, sample_name + "_" + group_id + ".png")
    os.makedirs(os.path.dirname(_output_plot), exist_ok=True)
    plt.savefig(_output_plot, dpi = 300)

    plt.close()
