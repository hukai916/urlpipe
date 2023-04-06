#!/usr/bin/env python

"""
To plot the UMI distribution pattern.

"""

import sys
import os
import pandas as pd

csv         = sys.argv[1]
umi_cutoff  = sys.argv[2]
outfile     = sys.argv[3]

if os.path.getsize(csv) > 0:
    col_names = ["umi", "read_count", "repeat_length"]
    df = pd.read_csv(csv, header = None, names = col_names)
    df_umi = df[df['read_count'] >= int(umi_cutoff)]
    df_res = df_umi["repeat_length"].value_counts().to_frame().reset_index()
    df_res = df_res.rename(columns = {'index': 'repeat_length', 'repeat_length':'count'})
    df_res = df_res.sort_values(by = "repeat_length")

    df_res.to_csv(outfile, index = False)
else:
    print("CSV file is empty!")
