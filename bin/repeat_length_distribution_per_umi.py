#!/usr/bin/env python

"""
To visualize the repeat length distribution by UMI group.

"""

import sys
import os
import random

csv = sys.argv[1]
cutoff      = sys.argv[2]
outfile_csv = sys.argv[3]
group_num   = sys.argv[4]

# read csv into dict
d_umi = {}
for line in open(csv):
    umi, length = line.strip().split(",")
    umi = umi.split("_")[1]
    if not umi in d_umi:
        d_umi[umi] = [length]
    else:
        d_umi[umi].append(length)

# shuffle d_umi
items = list(d_umi.items())
random.shuffle(items)
d_umi = dict(items)

d_cutoff = {}
for x in d_umi:
    if len(d_umi[x]) == int(cutoff):
        if not len(d_cutoff) > int(group_num):
            d_cutoff[x] = d_umi[x]
        else:
            break

# save to csv
with open(outfile_csv, "w") as f:
    for k in d_cutoff:
        for c in d_cutoff[k]:
            f.write(k + "," + str(c) + "\n")
