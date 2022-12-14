#!/usr/bin/env python

"""
To plot SD of read length distribution per UMI cutoffs.
Also plot the violin plot of read length distribution.

"""

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import os

sample_name = sys.argv[1]
outdir      = sys.argv[2]
cutoffs_str = sys.argv[3]
cutoff_0    = sys.argv[4]
cutoffs_file= sys.argv[5:]

# For STD scatter plot:
output_sd_plot = os.path.join(outdir, sample_name + "_std.png")
input_files = [cutoff_0] + cutoffs_file
sd_list = []
raw_list = []
_tem = []

for file in input_files:
    if os.path.getsize(file) > 0:
        tem = pd.read_csv(file, header = None)
        tem = tem[tem.iloc[:,0] != "problem"]
        tem = tem[tem.iloc[:, 0] != "plus"]
        data_tem = np.repeat(tem.iloc[:, 0], tem.iloc[:, 1])
        ls_tem = [int(x) for x in data_tem]
        raw_list.append(ls_tem)
        sd_list.append(np.std(ls_tem))
        if file.endswith("stat.csv"):
            # _tem.append("Cutoff_0")
            pass
        else:
            _cutoff = file.split(".csv")[0].split("_")[-1]
            _tem.append("Cutoff_" + _cutoff)
        # file: 5d_merge_repeat_dist_umi_correct/read_length_distribution/cutoff_100/mode/stat_mode_9_cutoff_100.csv


# X = ["No_correction", "Cutoff_1", "Cutoff_3", "Cutoff_10", "Cutoff_30", "Cutoff_100"]
# _tem = ["Cutoff_" + x.strip() for x in cutoffs_str.split(",")]
X = ["No_correction"] + _tem
Y = [round(x, 2) for x in sd_list]
print(X)
print(sd_list)
plt.scatter(X, sd_list)
plt.ylim(1, max(sd_list) + 2)
locs, labels = plt.xticks()
for i, label in enumerate(Y):
    plt.annotate(label, (X[i], Y[i] * 1.05))
fig = plt.gcf()
fig.patch.set_facecolor('xkcd:white')
plt.savefig(output_sd_plot)
plt.clf()

# For violin plot:
output_violin_raw_plot = os.path.join(outdir, sample_name + "_violin_raw.png")
output_violin_zoom_plot = os.path.join(outdir, sample_name + "_violin_zoom.png")
try:
    positions = [x + 1 for x in range(len(X))]
    plt.violinplot(raw_list, positions = positions)
    fig = plt.gcf()
    fig.patch.set_facecolor('xkcd:white')
    locs, labels = plt.xticks()
    plt.xticks(ticks = locs, labels = [""] + X + [""], rotation = 10)

    plt.ylim(0, 200)
    plt.savefig(output_violin_raw_plot)
    plt.ylim(130, 180)
    plt.savefig(output_violin_zoom_plot)
    plt.clf()
except:
    pass
