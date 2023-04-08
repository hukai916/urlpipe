#!/usr/bin/env python

"""
output: UMI,read count
"""

import sys
import pandas as pd

umi_reads  = sys.argv[1]
umi_cutoff = int(sys.argv[2])
outfile    = sys.argv[3]

df = pd.read_csv(umi_reads, header = None, names = ["umi", "read_count"])
df_res = df[df["read_count"] >= umi_cutoff]
res = df_res["read_count"].sum()

with open(outfile, "w") as f:
    f.write(str(res))