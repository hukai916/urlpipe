#!/usr/bin/env python
#%%
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
include_indel = sys.argv[6]
output_html = sys.argv[7]

# bins = ast.literal_eval("[(0,50), (51,60), (61,137), (138,154), (155,1000)]")
# use_ratio = "yes"
# use_repeat_unit_bp = "no"
# repeat_unit_bp = 3
# input_csv = "../test/b1b214983958d0cd2a4f519bdb43df/master_table_repeat_bin_umi_1.csv"
# include_indel = "yes"
# output_html = "../test/b1b214983958d0cd2a4f519bdb43df/test.html"
#%%

df = pd.read_csv(input_csv)
# %%
if use_repeat_unit_bp.lower() == "yes":
    if include_indel == "no":
        top_labels = [str(round(x[0]/repeat_unit_bp)) + "-" + str(round(x[1]/repeat_unit_bp)) for x in bins]
    else:
        top_labels = [str(round(x[0]/repeat_unit_bp)) + "-" + str(round(x[1]/repeat_unit_bp)) for x in bins if not x == "indel"]
        top_labels = ["indel"] + top_labels
else:
    if include_indel == "no":
        top_labels = [str(x[0]) + "-" + str(x[1]) for x in bins]
    else:
        top_labels = [str(x[0]) + "-" + str(x[1]) for x in bins if not x == "indel"]
        top_labels = ["indel"] + top_labels

if include_indel == "yes":
    bins = ["indel"] + bins

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
# y_data = df.columns.tolist()
# y_data = [x for x in y_data if not x == "repeat_length"]
y_data = df["sample_name"].tolist()

# %% 
x_data = []
for s in y_data:
    res = []
    for bin in bins:
        # accum_sum = df.loc[(df["repeat_length"] >= bin[0]) & (df["repeat_length"] < bin[1]), s].sum()
        df_tem = df[df["sample_name"] == s]
        if not bin == "indel":
            col_bin = str(bin[0]) + "-" + str(bin[1])
        else:
            col_bin = "indel"
        if include_indel == "no":
            df_tem.loc[:, "total"] = df_tem["total"] - df_tem["indel"]
        
        if use_ratio.lower() == "no":
            res.append(df_tem[col_bin].tolist()[0])
        else:
            res.append((df_tem[col_bin]/df_tem["total"]).tolist()[0])
    x_data.append(res)

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
                            x = 0, 
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

# Update xaxis range:
if use_ratio == "yes":
    xaxis = dict(range = [0, 1])
else:
    upper = max([item for sublist in x_data for item in sublist])
    xaxis = dict(range = [0, upper])

fig.update_layout(xaxis_title = stat,
                #   yaxis_title = "sample id",
                  title = "Repeat Length Distribution",
                  xaxis = xaxis)

fig.write_html(output_html)
# %%
