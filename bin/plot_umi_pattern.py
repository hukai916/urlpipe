#!/usr/bin/env python

"""
Plot UMI pattern: "reads per UMI" vs "log2(count)"
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

tsv = sys.argv[1]
sample_name = sys.argv[2]
outfile = sys.argv[3]

data=pd.read_csv(tsv, sep='\t')
x = data.iloc[:,0].to_frame()
weight = np.log2(data[["count"]])

# n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
n, bins, patches = plt.hist(x, weights = weight, bins=int(x.iloc[-1]))

plt.title(sample_name)
plt.xlabel("reads per UMI")
plt.ylabel("log2(count)")
plt.savefig(outfile, dpi = 300)
plt.close()
