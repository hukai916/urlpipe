#!/usr/bin/env python

"""
To plot the per mapped site read quality distribution.

Usage:
    python plot_per_site_qc.py per_site_*.csv outfile.html
"""

import sys
import os
import plotly.graph_objs as go
import plotly.offline as pyo
import numpy as np

per_site_qc_csv = sys.argv[1:-1]
out_html = sys.argv[-1]

qc_list = []
sn_list = [] # sample name

for file in per_site_qc_csv:
    with open(file, "r") as f:
        for line in f:
            tem = line.strip().split(",")
            if len(tem) > 1:
                qc_list.append([float(x) if x != "" else 0 for x in tem[1:]])
            else:
                qc_list.append([np.nan])
            basename = os.path.basename(file)
            sample_name = os.path.splitext(basename)[0].replace("mean_qc_", "")
            sn_list.append(sample_name)
            break # only one line

# Create the box plot
fig = go.Figure()
for i, lst in enumerate(qc_list):
    name = sn_list[i]
    fig.add_trace(go.Box(y = lst, name = name))

# Update the layout
fig.update_layout(
    title='Mean read quality distribution per category',
    xaxis = dict(title = 'Read category'),
    yaxis = dict(title = 'Read mean quality')
)

pyo.plot(fig, filename = out_html)
