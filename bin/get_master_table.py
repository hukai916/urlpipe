#!/usr/bin/env python

import sys
import pandas as pd

frac_csv = sys.argv[1:-3]
indel_csv = sys.argv[-3]
allele_number = int(sys.argv[-2])
outfile = sys.argv[-1]

# frac_csv = "repeat_length_fraction_umi_1.csv"
# indel_csv = "read_count_umi_cutoff_1.csv"
# allele_number = 1
# outfile = "test.csv"

df_res = pd.DataFrame()
for csv in frac_csv:
    df_tem = pd.read_csv(csv, dtype = {"sample_name": str})
    df_res = pd.concat([df_res, df_tem])

df_indel = pd.read_csv(indel_csv)

# add the indel column
indel = []
for sample in df_res["sample_name"]:
    indel.append(df_indel[str(sample)][0])
df_res["indel"] = indel

# ad the total and frac columns
if allele_number == 1:
    df_res["total"] = df_res["below"] + df_res["allele_1"] + df_res["above"] + df_res["indel"]
    df_res["below_frac"] = df_res["below"] / df_res["total"]
    df_res["allele_1_frac"] = df_res["allele_1"] / df_res["total"]
    df_res["above_frac"] = df_res["above"] / df_res["total"]
elif allele_number == 2:
    df_res["total"] = df_res["below"] + df_res["allele_1"] + df_res["between"] + df_res["allele_2"] + df_res["above"] + df_res["indel"]
    df_res["below_frac"] = df_res["below"] / df_res["total"]
    df_res["allele_1_frac"] = df_res["allele_1"] / df_res["total"]
    df_res["between_frac"] = df_res["between"] / df_res["total"]
    df_res["allele_2_frac"] = df_res["allele_2"] / df_res["total"]
    df_res["above_frac"] = df_res["above"] / df_res["total"]
df_res["indel_frac"] = df_res["indel"] / df_res["total"]

# reorder the colnames:
if allele_number == 1:
    new_column_order = ["sample_name", "start_allele_1", "end_allele_1", "total", "below", "allele_1", "above", "indel", "below_frac", "allele_1_frac", "above_frac", "indel_frac"]
elif allele_number == 2:
    new_column_order = ["sample_name", "start_allele_1", "end_allele_1", "start_allele_2", "end_allele_2", "total", "below", "allele_1", "between", "allele_2", "above", "indel", "below_frac", "allele_1_frac", "between_frac", "allele_2_frac", "above_frac", "indel_frac"]
df_res = df_res.reindex(columns = new_column_order)

with open(outfile, "w") as f:
    # f.write("test.csv")
    df_res.to_csv(outfile, index = False)