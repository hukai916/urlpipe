#!/usr/bin/env python

"""
Plot below/above fraction, mean, and std for all_sample.csv
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv = sys.argv[1]
outfile1 = sys.argv[2]
outfile2 = sys.argv[3]

data=pd.read_csv(csv, sep=',')

# Plot frac bar plot:
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
plt.savefig(outfile1, dpi = 300, bbox_inches = "tight")
plt.close()



# Plot repeat length bar plot:
N = len(labels)
ind = np.arange(N)  # the x locations for the groups
width = 0.4       # the width of the bars
below_mean = list(data.iloc[:,3])
below_std  = list(data.iloc[:,4])
above_mean = list(data.iloc[:,7])
above_std = list(data.iloc[:,8])

fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(ind, below_mean, width, yerr=below_std)
rects2 = ax.bar(ind+width, above_mean, width, yerr=above_std)

ax.set_ylabel('Repeat Length')
ax.set_title('Below and Above Cutoff Read Length by Sample')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(labels)

# ax.legend(('Below_frac', 'Above_frac') )
plt.xticks(rotation = 85)
plt.savefig(outfile2, dpi = 300, bbox_inches = "tight")
plt.close()
