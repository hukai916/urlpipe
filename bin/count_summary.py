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

df_res = df[["sample_name", "indel_count", "below_count", "between_count", "above_count"]]
df_res = df_res.rename(columns={'below_count': 'below',
                                'between_count': 'allele_1',
                                'above_count': 'above',
                                'indel_count': 'indel'
                                })

df_res['total'] = df_res.sum(axis=1)
for col in df_res.columns[:-1]:
    df_res[col + '_frac'] = df_res[col] / df_res['total']

cols_in_order = ["indel", "indel_frac", "below", "below_frac", "allele_1", "allele_1_frac", "above", "above_frac", "total"]
df_res = df_res[cols_in_order]
df_res.to_csv(outfile, index = False)
