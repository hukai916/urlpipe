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

df = pd.read_csv(csv, header = None)
df.columns = ["UMI", "Repeat_length"]
grouped_df = df.groupby('UMI')['Repeat_length']

fig = go.Figure()

for UMI, Repeat_length in grouped_df:
    fig.add_trace(go.Histogram(x = Repeat_length, name = UMI))

fig.update_layout(title='Repeat Length Distribution per Selected UMI group', xaxis_title = 'Repeat Length', yaxis_title = 'Count')

pyo.plot(fig, filename = outfile_html)
