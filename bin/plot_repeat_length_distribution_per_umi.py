#!/usr/bin/env python

"""
To visualize the repeat length distribution by UMI group.

"""

import sys
import os
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

csv = sys.argv[1]
cutoff      = sys.argv[2]
outfile_html = sys.argv[3]

df = pd.read_csv(csv, )
df_res.to_csv(outfile_csv, index = False)

# 2. output repeat_length_count_default_umi_0.html
    # create bar charts using frequencies
bars = [go.Bar(x = df["repeat_length"],
               y = df[col_sample],
               name = col_sample,
               width = 2,
               marker = dict(opacity = 0.5),
               visible = True if i == 0 else "legendonly"
               ) for i, col_sample in enumerate(df.columns[1:])]

fig = go.Figure(data = bars)

fig.update_layout(title = "Repeat Length Count Distribution",
                  xaxis_title = "Repeat Length (nt)",
                  yaxis_title = "Read Count")
outfile_html = "test.html"
# Save the plot as an interactive HTML file
pyo.plot(fig, filename = outfile_html)
