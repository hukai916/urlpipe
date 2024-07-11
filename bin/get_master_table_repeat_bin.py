#!/usr/bin/env python
#%%
import sys
import pandas as pd
import ast

repeat_length_count_csv = sys.argv[1] # repeat_length, sample1, sample2, ...
indel_csv = sys.argv[2] # sample1, sample2, ..., sampleX order should match repeat_length_count_csv
repeat_bins_str = sys.argv[3]
outfile = sys.argv[4]
#%% 
# repeat_length_count_csv = "../test/c9977df86539c34f71725affb71755/repeat_length_count_default_umi_3.csv"
# indel_csv = "../test/c9977df86539c34f71725affb71755/read_count_umi_cutoff_3.csv"
# repeat_bins_str = "[(0,50), (51,60), (61,137), (138,154), (155,1000)]"
# outfile = "../test/c9977df86539c34f71725affb71755/test.csv"

# read in repeat_length table: repeat_length, sample1, samplex, ...
df_res = pd.read_csv(repeat_length_count_csv)
# read in indel table: sample1, samplex, ...
df_indel = pd.read_csv(indel_csv)
# read in repeat bin ranges: [(0, 50), (xx, xx), ...]
repeat_bins = ast.literal_eval(repeat_bins_str)

#%%
# Sum up each repeat bin range for every sample:
def sum_within_range(df, col, lower, upper):
    return df.loc[(df["repeat_length"] >= lower) & (df["repeat_length"] <= upper), col].sum()

columns_to_sum = df_res.columns.difference(["repeat_length"])

#%%
summary_data = {}
# Calculate sums for each range and each column:
for lower, upper in repeat_bins:
    range_name = f"{lower}-{upper}"
    summary_data[range_name] = {}
    for col in columns_to_sum:
        summary_data[range_name][col] = sum_within_range(df_res, col, lower, upper)

# Create a new DataFrame from the summary data and add indel:
df_summary = pd.DataFrame(summary_data).rename_axis("sample_name")
df_indel_summary = df_indel.T.rename(columns={df_indel.T.columns[0]: "indel"}).rename_axis("sample_name")

df_final = pd.merge(df_summary, df_indel_summary, left_index = True, right_index = True, how = 'outer')

# Add total and percentage columns:
df_final['total'] = df_final.sum(axis = 1)
for col in df_final.columns:
    if col != 'total':
        df_final[f'{col}_percent'] = (df_final[col] / df_final['total']).round(3)
df_final.reset_index(inplace=True)

# Reorder columns
desired_columns = ["sample_name", "total"] + [col for col in df_final.columns if col not in ["sample_name", "total"]]
df_final = df_final[desired_columns]

with open(outfile, "w") as f:
    df_final.to_csv(outfile, index = False)