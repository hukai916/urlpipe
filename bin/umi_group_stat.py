#!/usr/bin/env python

"""
Summarize UMI group stat: UMI read_count mean mode ld

# ld: least distance determined corrected repeat length
"""

import sys
import os
from collections import Counter
from statistics import mean
import random
from collections import OrderedDict
from operator import itemgetter

csv = sys.argv[1]
sample_name = sys.argv[2]
outdir      = sys.argv[3]

# read csv into dict
d_umi = {}
for line in open(csv):
    umi, length = line.strip().split(",")
    umi = umi.split("_")[1]
    if str(length).isdigit():
        if not umi in d_umi:
            d_umi[umi] = [int(length)]
        else:
            d_umi[umi].append(int(length))

# save to stat
output_stat = os.path.join(outdir, sample_name + ".csv")
os.makedirs(os.path.dirname(output_stat), exist_ok=True)

with open(output_stat, "w") as f:
    for k in d_umi:
        ## Calculate the "mode" per UMI group
        mode = 0
        _data = Counter(d_umi[k])
        mode_freq = _data.most_common()[0][1]  # Returns the highest occurring item frequency
        mode_ls = []
        for x in _data.most_common(): # in terms of tie, randomly choose one
            if x[1] == mode_freq:
                mode_ls.append(x[0])
        mode = random.choice(mode_ls)

        ## Calculate the "ld" per UMI group: when tie, take mean of the ld_ls
        ld = 0
        ld_ls = []
        _tem = {} # store the height and its distance
        for height in _data: # height is the count, distance is nucleotide distance
            distance = 0
            for height2 in _data:
                distance += _data[height2] * abs(int(height2) - int(height)) # height * distance
            _tem[height] = distance
        sorted_x = OrderedDict(sorted(_tem.items(), key=itemgetter(1)))
        min_d    = list(sorted_x.values())[0]
        for height in sorted_x:
            if sorted_x[height] == min_d:
                ld_ls.append(height)
        if len(ld_ls) == 1:
            ld = ld_ls[0]
        else:
            ld = mean(ld_ls)

        # print(d_umi[k])
        res = ",".join([k, str(len(d_umi[k])), str(mean(d_umi[k])), str(mode), str(ld)]) + "\n"
        f.write(res)
