#!/usr/bin/env python

"""
Summarize UMI group stat: UMI read_count mean mode ld lsd

# ld: least distance determined corrected repeat length
# lsd: least squared distance which penalize larger distance even more
"""
#%%
import sys
from collections import Counter
from statistics import mean
import random
from collections import OrderedDict
from operator import itemgetter

csv = sys.argv[1]
umi_correction_method = sys.argv[2]
outfile = sys.argv[3]

#%%
def get_mode(d_umi_k):
    """
    When tie, randomly choose one mode.
    """
    _data = Counter(d_umi_k)
    mode_freq = _data.most_common()[0][1]  # Returns the highest occurring item frequency
    mode_ls = []
    for x in _data.most_common(): # in terms of tie, randomly choose one
        if x[1] == mode_freq:
            mode_ls.append(x[0])
    return(",".join([k, str(len(d_umi[k])), str(random.choice(mode_ls))]))

#%%
def get_least_distance(d_umi_k):
    """
    take the median, if tie, take the smaller value
    """
    ld_ls = []
    _tem = {} # store the height and its distance
    _data = Counter(d_umi_k) # repeat_length: count
    for height in _data: # height is the repeat_length, _data[height] is the count number
        distance = 0
        for height2 in _data:
            distance += _data[height2] * abs(int(height2) - int(height)) # count_number * distance
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
    return(",".join([k, str(len(d_umi[k])), str(ld)]))
#%%

def get_least_squared_distance(d_umi_k):
    """
    take the median, if tie, take the smaller value
    """
    ld_ls = []
    _tem = {} # store the height and its distance
    _data = Counter(d_umi_k) # repeat_length: count
    for height in _data: # height is the repeat_length, _data[height] is the count number
        distance = 0
        for height2 in _data:
            distance += _data[height2] * (int(height2) - int(height)) ** 2# count_number * distance
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
    return(",".join([k, str(len(d_umi[k])), str(ld)]))
#%%
def get_mean(d_umi_k):
    return(",".join([k, str(len(d_umi_k)), str(mean(d_umi_k))]))

#%%
# read csv into dict
d_umi = {} # stores umi as key and list of repeat lengths
for i, line in enumerate(open(csv)):
    if not i == 0:
        umi, length = line.strip().split(",")
        umi = umi.split("_")[1]
        if str(length).isdigit():
            if not umi in d_umi:
                d_umi[umi] = [int(length)]
            else:
                d_umi[umi].append(int(length))

with open(outfile, "w") as f:
    for k in d_umi:
        if umi_correction_method == "least_distance":
            res = get_least_distance(d_umi[k])
        elif umi_correction_method == "mode":
            res = get_mode(d_umi[k])
        elif umi_correction_method == "mean":
            res = get_mean(d_umi[k])
        elif umi_correction_method == "least_squared_distance":
            res = get_least_squared_distance(d_umi[k])
        f.write(res + "\n")
