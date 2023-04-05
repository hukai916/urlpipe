#!/usr/bin/env python

import sys
from utils import get_sample_name
import pandas as pd

csv = sys.argv[1:-1]
outfile_csv = sys.argv[-1]

sample_names = [get_sample_name(x, "repeat_length_count_", ".csv") for x in csv]

df_res = pd.DataFrame("repeat_length": [])

for i in range(len(csv)):
    df = pd.read_csv(csv[i])
    df = df.rename(columns={'count': sample_name[i]})
    df_res = pd.merge(df_res, df, on='repeat_length', how='outer')
    df_res.fillna(0, inplace = True)

df_res.to_csv(outfile_csv, index = False)
