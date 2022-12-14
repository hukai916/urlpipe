#!/usr/bin/env python

"""
Summarize the counts for reads: indel, below, between(unmodified), above

"""

import sys
import os
import pandas as pd

frac    = sys.argv[1]
indel   = sys.argv[2]
outfile = sys.argv[3]

df_frac  = pd.read_csv(frac)
df_indel = pd.read_csv(indel)

df = df_frac.join(df_indel, lsuffix='sample_name', rsuffix='sample_name')
df = df.loc[:,~df.columns.duplicated()]
df = df.rename(columns={"sample_namesample_name": "sample_name", "read_count": "indel_count"})

df_res = df[["sample_name", "indel_count", "blew_count", "between_count", "above_count"]]
df_res.to_csv(outfile, index = False)
