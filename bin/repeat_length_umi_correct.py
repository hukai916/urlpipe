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
umi_correction_method = sys.argv[2]
outfile = sys.argv[3]

def get_mode(d_umi):
    """
    When tie, randomly choose one mode.
    """
    mode = []
    for k in d_umi:
        _data = Counter(d_umi[k])
        mode_freq = _data.most_common()[0][1]  # Returns the highest occurring item frequency
        mode_ls = []
        for x in _data.most_common(): # in terms of tie, randomly choose one
            if x[1] == mode_freq:
                mode_ls.append(x[0])
        mode.append(",".join([k, str(len(d_umi[k])), str(random.choice(mode_ls))]))
    return(mode)

def get_least_distance(d_umi):
    """
    take the median, if tie, take the smaller value
    """
    res_ls = []
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
        # ld = mean(ld_ls) # mean is not a good choice since it creates many float values for downstream analysis
        ld = sorted(ld_ls)[int(len(ld_ls)/2)] # take the median, if tie, take the smaller value
    res_ld.append(",".join([k, str(len(d_umi[k])), str(ld)]))
    return(res_ld)

def get_mean(d_umi):
    res = []
    for k in d_umi:
        res.append(",".join([k, str(len(d_umi[k])), str(mean(d_umi[k]))]))

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

with open(outfile, "w") as f:
    if umi_correction_method == "least_distance":
        res = get_least_distance(d_umi)
    elif umi_correction_method == "mode":
        res = get_mode(d_umi)
    elif umi_correction_method == "mean":
        res = get_mean(d_umi)

    for line in res:
        f.write(line + "\n")
