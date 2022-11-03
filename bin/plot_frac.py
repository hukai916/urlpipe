#!/usr/bin/env python

"""
Plot below/above fraction, mean, and std for all_sample.csv
for 4d, the sample names are like: 13-PN.stat.csv
for 5d, the sample names are like: stat_mode_11-PN_cutoff_100.csv

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from pathlib import Path

csv = sys.argv[1]
outfile1 = sys.argv[2] # bar plot_plot
outfile2 = sys.argv[3] # violin plot
outfile2_raw = outfile2.replace(".png", ".raw.png")
outfile2_zoom = outfile2.replace(".png", ".zoom.png")
umi_cutoff   = sys.argv[4] # if 0, means raw without UMI correction

data=pd.read_csv(csv, sep=',')

fig, (ax1, ax2, ax3) = plt.subplots(1,3)
fig.set_size_inches(15,15)
fig.patch.set_facecolor('xkcd:white')

labels = list(data.iloc[:,0])
below_mean = list(data.iloc[:,2])
above_mean = list(data.iloc[:,6])
width = 0.35

# above and below
ax1.bar(labels, below_mean, width, label='below_frac')
ax1.bar(labels, above_mean, width, bottom=below_mean, label='above_frac')
ax1.set_ylabel('Fraction')
ax1.set_title('Below and Above Cutoff Reads Fraction by Sample')
ax1.legend()
ax1.set_xticklabels(labels, rotation = 85)

# above only
ax2.bar(labels, above_mean, width, label='above_frac')
ax2.set_ylabel('Fraction')
ax2.set_title('Above')
ax2.legend()
ax2.set_xticklabels(labels, rotation = 85)

# below only
ax3.bar(labels, below_mean, width, label='below_frac')
ax3.set_ylabel('Fraction')
ax3.set_title('Below')
ax3.legend()
ax3.set_xticklabels(labels, rotation = 85)

plt.savefig(outfile1, dpi = 600)
plt.clf()

# Violin plots
if umi_cutoff == 0:
    path = r'*.stat.csv'
    files = sorted(glob.glob(path))
    path = r'*cutoff_*.csv'
    files2 = sorted(glob.glob(path))
    files = files + files2
else:
    path = '*cutoff_' + str(umi_cutoff) + ".csv"
    #r'*cutoff_*.csv'
    path = path.encode('unicode_escape').decode() # make it raw string
    files = sorted(glob.glob(path))

raw_list = []
sample_list = []

for x in files:
    tem = pd.read_csv(x)
    tem = tem[tem.iloc[:,0] != "problem"]
    tem = tem[tem.iloc[:, 0] != "plus"]
    data_tem = np.repeat(tem.iloc[:, 0], tem.iloc[:, 1])
    ls_tem = [int(x) for x in data_tem]
    raw_list.append(ls_tem)

    p = Path(x)
    sample = p.stem.replace(".stat", "")
    sample_list.append(sample)

try:
    plt.violinplot(raw_list, positions = list(range(1, len(sample_list) + 1)),
                   vert = False,
                   showextrema = False)
    fig = plt.gcf()
    fig.patch.set_facecolor('xkcd:white')
    fig.set_size_inches(6, 6)

    plt.xlim(0, 200)
    locs, labels = plt.xticks()
    plt.yticks(ticks = list(range(1, len(sample_list) + 1)), labels = sample_list, fontsize = 8)
    plt.tick_params(axis="y",direction="in", pad=-130)
    plt.axvline(x = 151, color = 'r', linewidth = 1, linestyle = "dotted")
    plt.axvline(x = 166, color = 'r', linewidth = 1, linestyle = "dotted")
    plt.savefig(outfile2_raw, dpi = 600)
    plt.xlim(130, 180)
    plt.savefig(outfile2_zoom, dpi = 600)
except:
    pass
