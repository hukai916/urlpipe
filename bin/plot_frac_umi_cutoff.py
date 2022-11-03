#!/usr/bin/env python

"""
Plot UMI groups upon UMI correction using bargraph.

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from pathlib import Path

sample_name = sys.argv[1]
csv_no_correction = sys.argv[2]
umi_cutoff  = sys.argv[3]
outfile     = sys.argv[4]

umis = [x.strip() for x in umi_cutoff.split(",")]
res = [sum(pd.read_csv(csv_no_correction, header = None).iloc[:, 1])]
for umi in umis:
    file = "stat_mode_" + sample_name + "_cutoff_" + umi + ".csv"
    df = pd.read_csv(file, header = None)
    res.append(sum(df.iloc[:, 1]))
res[0] = int(res[0]/10)

fig, ax = plt.subplots()
width = 0.35
label = ["Cutoff_" + x for x in umis]
label = ["No_correction\n (count/10)"] + label

ax.bar(label, res, width)
ax.set_ylabel("Counts")
ax.set_title("UMI group count vs UMI cutoff")
plt.xticks(rotation = 10)

fig.savefig(outfile)
