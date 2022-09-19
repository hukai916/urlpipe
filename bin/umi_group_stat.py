#!/usr/bin/env python

"""
Summarize UMI group stat: UMI read_count mean mode
"""

import sys
import os
from collections import Counter
from statistics import mean

tsv = sys.argv[1]
sample_name = sys.argv[2]
outdir      = sys.argv[3]

# read tsv into dict
d_umi = {}
for line in open(tsv):
    umi, length = line.split()
    umi = umi.split("_")[1]
    if str(length).isdigit():
        if not umi in d_umi:
            d_umi[umi] = [length]
        else:
            d_umi[umi].append(length)

# save to stat
output_stat = os.path.join(outdir, sample_name + ".tsv")
os.makedirs(os.path.dirname(output_stat), exist_ok=True)
with open(output_stat, "w") as f:
    for k in d_umi:
        mode = 0
        _data = Counter(d_umi[k])
        mode_freq = _data.most_common()[0][1]  # Returns the highest occurring item frequency
        for x in _data.most_common():
            if x[1] == mode_freq:
                if mode > mode:
                    mode = x[0]
        print(d_umi[k])
        res = "\t".join([k, str(len(d_umi[k])), str(mean(d_umi[k])), str(mode)]) + "\n"
        f.write(res)
