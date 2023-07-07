import plotly.graph_objects as go
import numpy as np

def convert_to_phred(quality_string, offset = 33):
    phred_scores = []
    for char in quality_string:
        phred_score = ord(char) - offset
        phred_scores.append(phred_score)
    return phred_scores

# def add_scatter_trace(fig, x_values, y_values, trace_name = "test_trace"):
#     x = []
#     y = []
#     for i in range(len(x_values)):
#         x = x + [x_values[i]] * len(y_values[i])
#         y = y + list(y_values[i])
#     fig.add_trace(go.Scatter(x = x, y = y, mode = 'markers', marker = dict(size=3), name = trace_name))
#     return(fig)

test_csv = "per_site_qc_indel_5p_only_hs_snp1_BC10.csv"

sn = []
pos = []
qs = []

with open(test_csv, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        sn.append(tem[0])
        pos.append(int(tem[1]))
        qs.append(convert_to_phred(",".join(tem[2:])))



# x_values = []
# y_values = []
# for i, v in enumerate(pos):
#     x_values += [i] * len(qs[i])
#     y_values += qs[i]

# # Generate random data for demonstration
# x_values = np.linspace(0, 10, 10)
# y_values = np.random.randn(10, 35)  # 1000 x-values with 200 y-values each

# print(x_values)


# # Create a scatter plot for each x-value and y-values
fig = go.Figure()

for i, v in enumerate(pos):
    print(i)
    if i > 1000 and i < 1200:
        fig.add_trace(go.Box(x = v * np.ones_like(qs[i]), y = qs[i]))
    # if i > 100:
    #     break
    
# go.Box(x =np.ones_like(data_1), y=data_1, name='Group 1')


# print(x_values[1:100])
# print(y_values[1:100])

# fig.add_trace(go.Scatter(x = x_values[1:1000], y = y_values[1:1000], mode = 'markers', marker = dict(size = 3), name = "running"))
# add_scatter_trace(fig, x_values, y_values, "running")
# print(x_values)
    
fig.update_layout(
    title='Scatter Plot2',
    xaxis=dict(title='X'),
    yaxis=dict(title='Y'),
    showlegend = False
)

fig.show()