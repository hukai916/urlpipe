import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import os
import pandas as pd

def plot_repeat_dist(csv, output_file, sample_name, N, bin_number = "auto"):
    """
    csv:
        col1 col2
        repeat_length   count
    N = 5  # 1 tick every 5
    """
    # output plot:
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # We can set the number of bins with the *bins* keyword argument.
    df = pd.read_csv(csv, sep = ",", header = None)
    x, y = [], []
    _dict = {}

    for i in range(len(df.iloc[:,0])):
        _dict[str(df.iloc[i, 0])] = df.iloc[i, 1]


    bin_n = int(bin_number)
    # use bin_n or (int(df.iloc[-3, 0]) + 10)
    for i in list(range(bin_n - 2)) + ["plus", "problem"]:
        if not str(i) in _dict:
            x = x + [str(i)]
            y = y + [0]
        else:
            x = x + [str(i)]
            y = y + [_dict[str(i)]]

    # weight = [math.log2(x + 0.1) for x in df.iloc[:, 1]]
    weight = y
    # n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
    if bin_number == "auto":
        bin_n = len(x)
    else:
        bin_n = int(bin_number)
    n, bins, patches = plt.hist(x, weights = weight, bins=bin_n, range = (0, bin_n))

    # sparse the x-axis ticks:
    myticks = [i for i in range(bin_n) if not i%N]
    newlabels = [i for i in range(bin_n) if not i%N]

    if "plus" in list(x) and not "plus" in newlabels:
        # newlabels = newlabels[:-2]
        # myticks = myticks[:-2]
        newlabels = newlabels + ["+"]
        # myticks = myticks + [len(x) - 2]
        myticks = myticks + [bin_n + 1]

    plt.xticks(myticks, newlabels, rotation = 85)
    plt.tick_params(labelsize=8)
    plt.title(sample_name)

    plt.savefig(output_file, dpi = 300)
    plt.close()

def print_test():
    return("inside here")
