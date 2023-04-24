#!/usr/bin/env python

"""
To visualize the repeat length distribution by UMI group.

"""

import sys
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

csv = sys.argv[1]
cutoff      = sys.argv[2]
outfile_html = sys.argv[3]


# import os
# csv = "repeat_length_distribution_per_umi_1_07.csv"
# cutoff = 1
# outfile_html = "test.html"

try:
    df = pd.read_csv(csv, header = None)
except pd.errors.EmptyDataError:
    print("Error: CSV file (" + csv + ") is empty. Stopping.")
else:
    df.columns = ["UMI", "Repeat_length"]
    grouped_df = df.groupby('UMI')['Repeat_length']

    fig = go.Figure()

    x_range = [min(df["Repeat_length"]) - 1, max(df["Repeat_length"]) + 100]
    nbinsx  = x_range[1] - x_range[0] + 1
    x_ticks  = int(nbinsx / 5)

    for UMI, Repeat_length in grouped_df:
        fig.add_trace(go.Histogram(x = Repeat_length,
                                   name = UMI,
                                   nbinsx = nbinsx))

    fig.update_xaxes(
        range  = x_range,
        nticks = x_ticks
    )

    fig.update_layout(title = 'Repeat Length Distribution per Selected UMI group',           xaxis_title = 'Repeat Length',
                      yaxis_title = 'Count')

    pyo.plot(fig, filename = outfile_html)
