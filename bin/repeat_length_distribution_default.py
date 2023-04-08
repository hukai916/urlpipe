#!/usr/bin/env python

"""
Take repead_length_per_read_default_xxx.csv, generate diagnosis stat and output repeat_length_count_default_xxx.csv

"""

import pandas as pd
import sys

filename = sys.argv[1]
outfile_repeat_length_per_read = sys.argv[2]
outfile_repeat_length_count = sys.argv[3]
outfile_diagnosis_repeat_length_count = sys.argv[4]

df = pd.read_csv(filename, header = None, names=["read_id", "r1", "r2"])
# r1: repeat_length_by_r1

# R1 is number while R2 is nan:
df_count_r1_only = df[df["r1"].notnull() & df["r2"].isnull()]
count_r1_only = len(df_count_r1_only)

# R1 is nan while R2 is number:
df_count_r2_only = df[df["r1"].isnull() & df["r2"].notnull()]
count_r2_only = len(df_count_r2_only)

# boht R1 and R2 are numbers and R1 > R2:
df_count_r1_gt_r2 = df[df["r1"] > df["r2"]]
count_r1_gt_r2 = len(df_count_r1_gt_r2)

# both R1 and R2 are numbers and R1 == R2:
df_count_r1_eq_r2 = df[df["r1"] == df["r2"]]
count_r1_eq_r2 = len(df_count_r1_eq_r2)

# both R1 and R2 are numbers and R1 < R2:
df_count_r1_lt_r2 = df[df["r1"] < df["r2"]]
count_r1_lt_r2 = len(df_count_r1_lt_r2)

# both R1 and R2 are nan:
df_count_both_na = df[df["r1"].isnull() & df["r2"].isnull()]
count_both_na = len(df_count_both_na)

def df_rename(df):
    df = df.rename(columns = {"r1": "repeat_length", "r2": "repeat_length"})
    return(df)

df_r1_only = df_rename(df_count_r1_only[["read_id", "r1"]])
df_r2_only = df_rename(df_count_r2_only[["read_id", "r2"]])
df_r1_gt_r2 = df_rename(df_count_r1_gt_r2[["read_id", "r2"]])
df_r1_lt_r2 = df_rename(df_count_r1_lt_r2[["read_id", "r1"]])
df_r1_eq_r2 = df_rename(df_count_r1_eq_r2[["read_id", "r1"]])

df_res = pd.concat([df_r1_only, df_r2_only, df_r1_gt_r2, df_r1_lt_r2, df_r1_eq_r2])
df_res['repeat_length'] = df_res['repeat_length'].astype(int)

# output repeat length per read:
df_res.to_csv(outfile_repeat_length_per_read, index = False)

# output diagnosis:
with open(outfile_diagnosis_repeat_length_count, "w") as f:
    f.write("R1_length_R2_nan: " + str(count_r1_only) + "\n")
    f.write("R1_nan_R1_length: " + str(count_r2_only) + "\n")
    f.write("R1_eq_R2: " + str(count_r1_eq_r2) + "\n")
    f.write("R1_gt_R2: " + str(count_r1_gt_r2) + "\n")
    f.write("R1_lt_R2: " + str(count_r1_lt_r2) + "\n")
    f.write("R1_nan_R2_nan: " + str(count_both_na) + "\n")

# output repeat length count
counts = df_res['repeat_length'].value_counts()

counts_df = counts.to_frame()
counts_df = counts_df.rename(columns = {'repeat_length': 'count'})
counts_df = counts_df.sort_index().reset_index()
counts_df = counts_df.rename(columns={'index': 'repeat_length'})
counts_df.to_csv(outfile_repeat_length_count, index = False)
