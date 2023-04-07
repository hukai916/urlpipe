#!/usr/bin/env python

import sys
from utils import get_sample_name
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo


csv = sys.argv[1:-5]
outfile_csv = sys.argv[-4]
outfile_html = sys.argv[-3]
prefix = sys.argv[-2]
suffix = sys.argv[-1]

sample_names = [get_sample_name(x, prefix, suffix) for x in csv]

df_res = pd.DataFrame(columns = ["repeat_length"])

for i in range(len(csv)):
    df = pd.read_csv(csv[i])
    df = df.rename(columns={'count': sample_names[i]})
    df_res = pd.merge(df_res, df, on='repeat_length', how='outer')
    df_res.fillna(0, inplace = True)

df_res = df_res.sort_values('repeat_length')

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

fig.update_layout(title = "Repeat Length Count Distribution",
                  xaxis_title = "Repeat Length (nt)",
                  yaxis_title = "Read Count")

# Save the plot as an interactive HTML file
pyo.plot(fig, filename = outfile_html)
