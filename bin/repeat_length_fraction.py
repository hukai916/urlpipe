#!/usr/bin/env python

"""
Dev notes:
1.  pd.cut follows: (] range by default

"""

import sys
import pandas as pd
# from utils import get_std

allele_number  = int(sys.argv[1])
sample_name    = sys.argv[2]
repeat_length_count_table = sys.argv[3]
outfile        = sys.argv[4]
start_allele_1 = sys.argv[5]
end_allele_1   = sys.argv[6]
start_allele_2 = sys.argv[7]
end_allele_2   = sys.argv[8]

# allele_number = 2
# sample_name = "hs_BC22"
# repeat_length_count_table = "repeat_length_count_default_umi_30.csv"
# outfile = "test.csv"
# start_allele_1 = 50
# end_allele_1 = 60
# start_allele_2 = 140
# end_allele_2 = 155
# print(outfile)

# 1. read in sample column into df and add bin:
df = pd.read_csv(repeat_length_count_table)[["repeat_length", sample_name]]
empty = 0
if df.empty or df[sample_name].sum() == 0:
    empty = 1
    _total_count = 0
else:
    if allele_number == 1:        
        bins = [min(list(df["repeat_length"]) + [int(start_allele_1)]) - 1, int(start_allele_1), int(end_allele_1), max(list(df["repeat_length"]) + [int(end_allele_1)]) + 1]
    elif allele_number == 2:
        # print(df["repeat_length"] + int(start_allele_1))
        bins = [min(list(df["repeat_length"]) + [int(start_allele_1)]) - 1, int(start_allele_1), int(end_allele_1), int(start_allele_2), int(end_allele_2), max(list(df["repeat_length"]) + [int(end_allele_2)]) + 1]
    df["bin"] = pd.cut(df["repeat_length"], bins = bins)

 # 2. calculate count, fraction, mean, std:
    # first deal with edge case:
if allele_number == 1:
    header = "sample_name,start_allele_1,end_allele_1,below,allele_1,above"
elif allele_number == 2:
    header = "sample_name,start_allele_1,end_allele_1,start_allele_2,end_allele_2,below,allele_1,between,allele_2,above"
    _total_count = df[sample_name].sum()


if _total_count == 0 or empty:
    if allele_number == 1:
        res = header + "\n" + sample_name + ",0,0,0,0,0"
    elif allele_number == 2:
        res = header + "\n" + sample_name + ",0,0,0,0,0,0,0,0,0"
        print(res)
    with open(outfile, "w") as f:
        f.write(res + "\n")
else:
    # sys.exit("Total count is 0!")

    df["weighted_length"] = df["repeat_length"] * df[sample_name]
    _count = df.groupby("bin")[sample_name].sum()
    _fraction = _count / _total_count
    _mean = df.groupby("bin")["weighted_length"].sum() / df.groupby("bin")[sample_name].sum()
    # _std = df.groupby("bin")["repeat_length"].agg(lambda x: get_std(x))

    # 3. output to csv:
    if allele_number == 1:
        res = sample_name + "," + str(start_allele_1) + "," + str(end_allele_1) + ","
    elif allele_number == 2:
        res = sample_name + "," + str(start_allele_1) + "," + str(end_allele_1) + "," + str(start_allele_2) + "," + str(end_allele_2) + ","
        
    for i in range(len(_count)):
        # res = res + ",".join([str(_count.iloc[i]), str(_fraction.iloc[i]), str(_mean.iloc[i]), str(_std.iloc[i])])
        res = res + str(_count.iloc[i])
        if not i == len(_count) - 1:
            res = res + ","

    with open(outfile, "w") as f:
        f.write(header + "\n")
        f.write(res)
