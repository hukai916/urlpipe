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
cutoff_0    = sys.argv[3] # tem read length without UMI correction.
cutoff_1    = sys.argv[4]
cutoff_3    = sys.argv[5]
cutoff_10   = sys.argv[6]
cutoff_30   = sys.argv[7]
cutoff_100  = sys.argv[8]

# For STD scatter plot:
output_sd_plot = os.path.join(outdir, "plot_std_", sample_name + ".png")
input_files = [cutoff_0, cutoff_1, cutoff_3, cutoff_10, cutoff_30, cutoff_100]
sd_list = []
raw_list = []

for file in input_files:
    print(file)
    tem = pd.read_csv(file)
    tem = tem[tem.iloc[:,0] != "problem"]
    tem = tem[tem.iloc[:, 0] != "plus"]
    data_tem = np.repeat(tem.iloc[:, 0], tem.iloc[:, 1])
    ls_tem = [int(x) for x in data_tem]
    raw_list.append(ls_tem)
    sd_list.append(np.std(ls_tem))

X = ["No_correction", "Cutoff_1", "Cutoff_3", "Cutoff_10", "Cutoff_30", "Cutoff_100"]
Y = [round(x, 2) for x in sd_list]
plt.scatter(X, sd_list)
plt.ylim(1, 10)
locs, labels = plt.xticks()
for i, label in enumerate(Y):
    plt.annotate(label, (X[i], Y[i] * 1.05))
fig = plt.gcf()
fig.patch.set_facecolor('xkcd:white')
plt.savefig(output_sd_plot)

# For violin plot:
output_violin_raw_plot = os.path.join(outdir, "plot_violin_raw_", sample_name + ".png")
output_violin_zoom_plot = os.path.join(outdir, "plot_violin_zoom_", sample_name + ".png")

plt.violinplot(raw_list)
fig = plt.gcf()
fig.patch.set_facecolor('xkcd:white')
locs, labels = plt.xticks()
plt.xticks(ticks = locs, labels = [""] + X + [""], rotation = 10)

plt.ylim(0, 200)
plt.savefig(output_violin_raw_plot)
plt.ylim(130, 180)
plt.savefig(output_violin_zoom_plot)
