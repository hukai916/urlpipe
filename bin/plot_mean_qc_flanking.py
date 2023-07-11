#!/usr/bin/env python

"""
To plot the read mean quality flanking distribution.

Usage:
    python plot_mean_qc_flanking.py left_flanking_score.txt right_flanking_score.txt outfile.html outstat.csv sample_name read_id.csv
"""

import sys
import plotly.graph_objs as go
import plotly.offline as pyo
import numpy as np

left_flanking_qc = sys.argv[1]
right_flanking_qc = sys.argv[2]
outhtml = sys.argv[3]
outcsv  = sys.argv[4]
sample_name = sys.argv[5]
out_read_id_mean_qc = sys.argv[6]

left_c, right_c = 0, 0
left_na, right_na = 0, 0
left_mean, right_mean = [], []
read_id_mean_qc = {}

with open(left_flanking_qc, "r") as f:
    for line in f: 
        tem = line.strip().split(",")
        if not tem == "":
            left_c += 1
            if tem[1] == "NA" or tem[1] == "":
                left_na += 1
            else:
                _left_mean = np.mean([int(x) for x in tem[1:]])
                left_mean.append(_left_mean)
                read_id_mean_qc[tem[0]] = [_left_mean]
            
with open(right_flanking_qc, "r") as f:
    for line in f: 
        tem = line.strip().split(",")
        if not tem == "":
            right_c += 1
            if tem[1] == "NA" or tem[1] == "":
                right_na += 1
            else:
                _right_mean = np.mean([int(x) for x in tem[2:]])
                right_mean.append(_right_mean)
                if tem[0] in read_id_mean_qc:
                    read_id_mean_qc[tem[0]].append(_right_mean)

with open(out_read_id_mean_qc, "w") as f:
    for id in read_id_mean_qc:
        if len(read_id_mean_qc[id]) == 2:
            f.write(str(id) + "," + str(read_id_mean_qc[id][0]) + "," + str(read_id_mean_qc[id][1]) + "\n")

with open(outcsv, "w") as f:
    res = [sample_name, str(left_c), str(right_c), str(left_na), str(right_na), str(left_c / (left_c + left_na)), str(right_c / (right_c + right_na))]
    f.write(",".join(res) + "\n")

# Create the box plot
fig = go.Figure()
fig.add_trace(go.Violin(y = left_mean, name = "left_flanking", box_visible = True, meanline_visible = True))
fig.add_trace(go.Violin(y = right_mean, name = "right_flanking", box_visible = True, meanline_visible = True))

# Update the layout
fig.update_layout(
    title='Mean read quality for repeat flanking regions',
    xaxis = dict(title = 'Region'),
    yaxis = dict(title = 'Mean read quality')
)

pyo.plot(fig, filename = outhtml)
