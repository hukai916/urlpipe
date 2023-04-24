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

# prefix = ""
# suffix = ".csv"
# csv  = ['01.csv', '05.csv', '09.csv', '13.csv', '17.csv', '21.csv', '25.csv', '29.csv', '33.csv', '37.csv', '41.csv', '02.csv', '06.csv', '10.csv', '14.csv', '18.csv', '22.csv', '26.csv', '30.csv', '34.csv', '38.csv', '42.csv', '03.csv', '07.csv', '11.csv', '15.csv', '19.csv', '23.csv', '27.csv', '31.csv', '35.csv', '39.csv', '04.csv', '08.csv', '12.csv', '16.csv', '20.csv', '24.csv', '28.csv', '32.csv', '36.csv', '40.csv']
# col_name = "sample_id,read_count_merge,read_count_non_merge".split(",")
# outfile_csv = "test.csv"
# outfile_html = "test.html"

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