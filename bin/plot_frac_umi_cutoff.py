#!/usr/bin/env python

"""
Plot UMI groups upon UMI correction using bargraph.

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import glob

sample_name = sys.argv[1]
csv_no_correction = sys.argv[2]
umi_cutoff  = sys.argv[3]
outfile     = sys.argv[4]

umis = [x.strip() for x in umi_cutoff.split(",")]
res  = pd.read_csv(csv_no_correction, header = None).iloc[:, [2,6]] # above/below percentage columns
for x in umis:
    filename = sample_name + "_frac_*_cutoff_" + x + ".csv"
    file = glob.glob(filename)[0]
    res = pd.concat([res, pd.read_csv(file, header = None).iloc[:, [2,6]]], ignore_index = True)

fig, (ax1, ax2) = plt.subplots(1,2)
width = 0.35
label = ["Cutoff_" + x for x in umis]
label = ["No_correction"] + label

print(res)
print(csv_no_correction)

ax1.bar(label, res.iloc[:, 0], width)
ax1.set_ylabel("Percentage")
ax1.set_title("Above cutoff fraction")
ax1.set_xticklabels(label, rotation = 25, fontsize = 8)

ax2.bar(label, res.iloc[:, 1], width)
ax2.set_ylabel("Percentage")
ax2.set_title("Below fraction cutoff")
ax2.set_xticklabels(label, rotation = 25, fontsize = 8)

fig.tight_layout()

fig.savefig(outfile, dpi = 600)
