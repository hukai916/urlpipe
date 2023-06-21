#!/usr/bin/env python

import sys
from utils import get_sample_name
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

csv = sys.argv[1:-5]
outfile_csv = sys.argv[-5]
outfile_html = sys.argv[-4]
prefix = sys.argv[-3]
suffix = sys.argv[-2]
xlim   = sys.argv[-1]

# import os 
# # Get a list of all CSV files in the directory
# csv = [f for f in os.listdir("./") if f.endswith('.csv')]
# outfile_csv = "test.csv"
# outfile_html = "test.html"
# prefix = "repeat_length_count_"
# suffix = ".csv"

prefix = "repeat_length_count_default_"
suffix = "_umi_100.csv"
outfile_csv = "test.csv"
xlim = 500
outfile_html = "test.html"
csv = ['repeat_length_count_default_hs_BC01_umi_100.csv', 'repeat_length_count_default_hs_BC13_umi_100.csv', 'repeat_length_count_default_hs_BC02_umi_100.csv', 'repeat_length_count_default_hs_BC14_umi_100.csv', 'repeat_length_count_default_hs_BC03_umi_100.csv', 'repeat_length_count_default_hs_BC15_umi_100.csv', 'repeat_length_count_default_hs_BC04_umi_100.csv', 'repeat_length_count_default_hs_BC16_umi_100.csv', 'repeat_length_count_default_hs_BC05_umi_100.csv', 'repeat_length_count_default_hs_BC17_umi_100.csv', 'repeat_length_count_default_hs_BC06_umi_100.csv', 'repeat_length_count_default_hs_BC18_umi_100.csv', 'repeat_length_count_default_hs_BC07_umi_100.csv', 'repeat_length_count_default_hs_BC19_umi_100.csv', 'repeat_length_count_default_hs_BC08_umi_100.csv', 'repeat_length_count_default_hs_BC20_umi_100.csv', 'repeat_length_count_default_hs_BC09_umi_100.csv', 'repeat_length_count_default_hs_BC21_umi_100.csv', 'repeat_length_count_default_hs_BC10_umi_100.csv', 'repeat_length_count_default_hs_BC22_umi_100.csv', 'repeat_length_count_default_hs_BC11_umi_100.csv', 'repeat_length_count_default_hs_BC23_umi_100.csv', 'repeat_length_count_default_hs_BC12_umi_100.csv']

sample_names = [get_sample_name(x, prefix, suffix) for x in csv]

df_res = pd.DataFrame(columns = ["repeat_length"])

for i in range(len(csv)):
    df = pd.read_csv(csv[i])
    df = df.rename(columns={'count': sample_names[i]})
    df_res = pd.merge(df_res, df, on='repeat_length', how='outer')
    df_res.fillna(0, inplace = True)

df_res = df_res.sort_values('repeat_length')

# move 'repeat_length' as col0:
df_res.insert(0, 'repeat_length', df_res.pop('repeat_length'))

# 1. output repeat_length_count_default_umi_0.csv
df_res.to_csv(outfile_csv, index = False)

# 2. output repeat_length_count_default_umi_0.html
    # create bar charts using frequencies
bars = [go.Bar(x = df_res["repeat_length"],
               y = df_res[col_sample],
               name = col_sample,
               width = 2,
               marker = dict(opacity = 0.5),
               visible = True if i == 0 else "legendonly"
               ) for i, col_sample in enumerate(df_res.columns[1:])]

fig = go.Figure(data = bars)

# x_range = [min(df["repeat_length"]) - 1, min(int(xlim), max(df["repeat_length"]) + 100 )]
if df_res.empty:
    x_range = [0, min(int(xlim), 100 )]
else:
    x_range = [0, min(int(xlim), max(df["repeat_length"]) + 100 )]
    
nbinsx  = x_range[1] - x_range[0] + 1
x_ticks  = int(nbinsx / 5)

fig.update_xaxes(
    range  = x_range,
    nticks = x_ticks
)

fig.update_layout(title = "Repeat Length Count Distribution",
                  xaxis_title = "Repeat Length (nt)",
                  yaxis_title = "Read Count")

# Save the plot as an interactive HTML file
pyo.plot(fig, filename = outfile_html)
