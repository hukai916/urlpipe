#!/usr/bin/env python

import sys
from utils import get_sample_name
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

csv = sys.argv[1:-4]
outfile_csv = sys.argv[-4]
outfile_html = sys.argv[-3]
prefix = sys.argv[-2]
suffix = sys.argv[-1]

sample_names = [get_sample_name(x, prefix, suffix) for x in csv]

df_res = pd.DataFrame()

for i in range(len(csv)):
    df = pd.read_csv(csv[i], header = None, names = [sample_names[i]])
    df_res[sample_names[i]] = df[sample_names[i]]

# 1. output read_count_umi_cutoff_x.csv
df_res.to_csv(outfile_csv, index = False)

# 2. output read_count_umi_cutoff_x.html
    # create bar charts using frequencies
bars = [go.Bar(x = list(df_res.columns),
               y = list(df_res.iloc[0]),
               marker = dict(opacity = 0.5)
               )]

fig = go.Figure(data = bars)

fig.update_layout(title = "Read Count at UMI cutoff" + str(),
                  xaxis_title = "Sample Name",
                  yaxis_title = "Read Count")

# Save the plot as an interactive HTML file
pyo.plot(fig, filename = outfile_html)
