#!/usr/bin/env python

import sys
from utils import get_sample_name
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

csv = sys.argv[1:-5]
col_name = sys.argv[-5].split(",")
outfile_csv = sys.argv[-4]
outfile_html = sys.argv[-3]
prefix = sys.argv[-2]
suffix = sys.argv[-1]

sample_names = [get_sample_name(x, prefix, suffix) for x in csv]

df_res = pd.DataFrame()

for i in range(len(csv)):
    df = pd.read_csv(csv[i], header = None, dtype = {0: str}, names = col_name)
    df_res = pd.concat([df_res, df], ignore_index = True) 

df_res["fraction"] = df_res[col_name[1]] / (df_res[col_name[1]] + df_res[col_name[2]])

# 1. output *.csv
df_res.to_csv(outfile_csv, index = False)

# 2. output *.html
    # create bar charts using values from col "fraction"
bars = [go.Bar(x = list(df_res[col_name[0]]),
               y = list(df_res["fraction"]),
               marker = dict(opacity = 0.5)
               )]

fig = go.Figure(data = bars)

fig.update_layout(title = "Fraction for each sample",
                  xaxis_title = "Sample Name",
                  yaxis_title = "Fraction")

# Save the plot as an interactive HTML file
pyo.plot(fig, filename = outfile_html)