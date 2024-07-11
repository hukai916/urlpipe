#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import ast


bins = sys.argv[1]
bins = ast.literal_eval(bins)
use_ratio = sys.argv[2]
use_repeat_unit_bp = sys.argv[3]
repeat_unit_bp = int(sys.argv[4])
input_csv = sys.argv[5]
output_html = sys.argv[6]

df = pd.read_csv(input_csv)
# %%
if use_repeat_unit_bp.lower() == "yes":
    top_labels = [str(round(x[0]/repeat_unit_bp)) + "-" + str(round(x[1]/repeat_unit_bp)) for x in bins]
else:
    top_labels = [str(x[0]) + "-" + str(x[1]) for x in bins]

# %%
# generate color
def generate_colors(num_colors):
    color_palette = plt.cm.tab10.colors
    num_colors = min(num_colors, len(color_palette))
    selected_colors = color_palette[:num_colors]
    rgba_colors = [tuple(int(x * 255) for x in color[:3]) + (1.0,) for color in selected_colors]
    rgba_strings = [f'rgba{color}' for color in rgba_colors]
    return rgba_strings
    
colors = generate_colors(len(bins))

# %%
y_data = df.columns.tolist()
y_data = [x for x in y_data if not x == "repeat_length"]

# %% 
x_data = []
for s in y_data:
    res = []
    for bin in bins:
        accum_sum = df.loc[(df["repeat_length"] >= bin[0]) & (df["repeat_length"] < bin[1]), s].sum()
        res.append(accum_sum)
    x_data.append(res)
    
# %%
if use_ratio.lower() == "yes": # convert to ratio if specified
    # x_data_ratio = x_data.copy() # direct assign create a link that would be updated when updating x_data_ratio
    x_data_ratio = []
    for x in x_data:
        total_count = sum(x)
        x_data_ratio.append([round(c/sum(x) * 100, 1) for c in x])
    x_data = x_data_ratio

# %%
fig = go.Figure() 
x_data_t = np.array(x_data).transpose().tolist()
for i in range(0, len(x_data_t)):
    fig.add_trace(go.Bar(
        x=x_data_t[i], y=y_data,
        name = top_labels[i],
        orientation = 'h',
        marker=dict(
            color=colors[i],
            line=dict(color='rgb(248, 248, 249)', width=1)
        )
    ))

fig.update_layout(
    xaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=True,
        zeroline=False,
        domain=[0.15, 1]
    ),
    yaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=False,
        zeroline=False,
    ),
    barmode='stack',
    paper_bgcolor='rgb(248, 248, 255)',
    plot_bgcolor='rgb(248, 248, 255)',
    margin=dict(l=120, r=10, t=140, b=80),
    showlegend=True,
)

# Label the y-axis
annotations = []
for i, (yd, xd) in enumerate(zip(y_data, x_data)):
    annotations.append(dict(xref = 'x', 
                            yref = 'y',
                            x = 0.14, 
                            y = i,
                            xanchor = 'right',
                            text = str(yd),
                            showarrow = False, 
                            align = 'right'))

fig.update_layout(annotations=annotations)

# Update x, y-axis title, and plot title
if use_ratio.lower() == "yes":
    stat = "Ratio"
elif use_ratio.lower() == "no":
    stat = "Read count"
if use_repeat_unit_bp.lower() == "yes":
    unit = "repeat unit"
elif use_repeat_unit_bp.lower() == "no":
    unit = "base pair"
    
fig.update_layout(xaxis_title = stat,
                  yaxis_title = "sample id",
                  title = "Repeat Length Distribution")

fig.write_html(output_html)