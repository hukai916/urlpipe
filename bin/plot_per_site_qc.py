#!/usr/bin/env python

"""
To plot the per mapped site read quality distribution.

Usage:
    python plot_per_site_qc.py per_site_*.csv outfile.html
"""

import plotly.graph_objects as go
import plotly.offline as pyo
import numpy as np
import sys
import os

per_site_csv = sys.argv[1:-1]
per_site_range_to_display = sys.argv[-1]
range_start, range_end = [int(x) for x in per_site_range_to_display.split(":")]

def convert_to_phred(quality_string, offset = 33):
    phred_scores = []
    for char in quality_string:
        phred_score = ord(char) - offset
        phred_scores.append(phred_score)
    return phred_scores

for csv in per_site_csv:
    sn = []
    pos = []
    qs = []
    
    basename = os.path.basename(csv)
    sample_name = os.path.splitext(basename)[0].replace("per_site_qc_", "")

    with open(csv, "r") as f:
        for line in f:
            tem = line.strip().split(",")
            sn.append(tem[0])
            pos.append(int(tem[1]))
            qs.append(convert_to_phred(",".join(tem[2:])))

    fig = go.Figure()

    for i, v in enumerate(pos):
        if i >= range_start and i <= range_end:
            fig.add_trace(go.Box(x = v * np.ones_like(qs[i]), y = qs[i]))

    fig.update_layout(
        title='Read quality per mapped site',
        xaxis=dict(title = 'Coordinate in reference'),
        yaxis=dict(title = 'Per base quality'),
        showlegend = False
    )

    pyo.plot(fig, filename = sample_name + ".html")