#!/usr/bin/env python

"""
Plot below/above fraction, mean, and std for all_sample.csv
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv = sys.argv[1]
outfile = sys.argv[2]

data=pd.read_csv(csv, sep=',')

fig, ax = plt.subplots()
labels = list(data.iloc[:,0])
below_mean = list(data.iloc[:,2])
above_mean = list(data.iloc[:,6])
width = 0.35

ax.bar(labels, below_mean, width, label='below_frac')
ax.bar(labels, above_mean, width, bottom=below_mean,
       label='above_frac')
ax.set_ylabel('Fraction')
ax.set_title('Below and Above Cutoff Reads Fraction by Sample')
ax.legend()
plt.xticks(rotation = 85)
plt.savefig(outfile, dpi = 300, bbox_inches = "tight")
plt.close()
# x = data.iloc[:,0].to_frame()
# weight = np.log2(data[["count"]])
#
# # n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
# if bin_number == "auto":
#     bin_n = int(x.iloc[-1])
# else:
#     bin_n = int(bin_number)
# n, bins, patches = plt.hist(x, weights = weight, bins = bin_n, range = (0, bin_n))
#
# plt.title(sample_name)
# plt.xlabel("reads per UMI")
# plt.ylabel("log2(count)")
# plt.savefig(outfile, dpi = 300)
# plt.close()
