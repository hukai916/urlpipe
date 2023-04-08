#!/usr/bin/env python

"""
output: UMI,read count
"""

import sys
import pandas as pd

umi_reads  = sys.argv[1]
umi_cutoff = int(sys.argv[2])
outfile    = sys.argv[3]

df = pd.read_csv(umi_reads, header = None)
res = df[df.iloc[1] >= umi_cutoff]

with open(outfile, "w") as f:
    f.write(res)