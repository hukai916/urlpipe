import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import sys
import pandas as pd

samples = ['1-Striatum-9285-D10A-CAGr-5pR3', '2-Striatum-9285-NI', '3-Striatum-9289-D10A-CAGr-5pR3',
           '4-Striatum-9289-NI', '5-Striatum-9290-D10A-CAGr-5pR3', '6-Striatum-9290-NI',
           '7-Striatum-9291-NI', '8-Striatum-9291-D10A-CAGr-5pR3', '9-PN-D10A-CAGr-5pR3-R1',
           '10-PN-D10A-CAGr-5pR3-R2', '11-PN-D10A-CAGr-5pR3-R3', '12-PN-D10A-CAGr-5pR3-R4',
           '13-PN-D10A-CAGr-R1', '14-PN-D10A-CAGr-R2', '15-PN-D10A-CAGr-R3',
           '16-PN-D10A-CAGr-R4', '17-PN-D10A-5pR3-R1', '18-PN-D10A-5pR3-R2',
           '19-PN-D10A-5pR3-R3', '20-PN-D10A-5pR3-R4', '21-PN-D10A-Rosa-R1',
           '22-PN-D10A-Rosa-R2','23-PN-D10A-Rosa-R3', '24-PN-D10A-Rosa-R4']

colNum = 2
# fig, axs = plt.subplots(int(24 / colNum), colNum, sharey=True, tight_layout=True)
# fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

# We can set the number of bins with the *bins* keyword argument.
for k, sample in enumerate(samples):
    sample = samples[k]
    _i = k // colNum
    _j = k - colNum * (k // colNum)
    df = pd.read_csv("../res/count_cag/" + sample + "_R1_001.fastq.tsv", sep = "\t")
    x, y = [], []
    _dict = {}

    for i in range(len(df.iloc[:,0])):
        _dict[df.iloc[i, 0]] = df.iloc[i, 1]

    for i in list(range(int(df.iloc[-3, 0]) + 10)) + ["plus", "problem"]:
        if not str(i) in _dict:
            x = x + [str(i)]
            y = y + [0]
        else:
            x = x + [str(i)]
            y = y + [_dict[str(i)]]

    # weight = [math.log2(x + 0.1) for x in df.iloc[:, 1]]
    weight = y
    # n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
    n, bins, patches = plt.hist(x, weights = weight, bins=len(x), edgecolor='black')

    # sparse the x-axis ticks:
    N = 5  # 1 tick every 5
    # xticks_pos, xticks_labels = plt.xticks()  # get all axis ticks

    # myticks = [j for i,j in enumerate(xticks_pos) if not i%N]  # index of selected ticks
    # newlabels = [i for i,label in enumerate(xticks_labels) if not i%N]
    myticks = [i for i in range(len(x)) if not i%N]
    newlabels = [i for i in range(len(x)) if not i%N]

    if "plus" in list(x) and not "plus" in newlabels:
        # newlabels = newlabels + ["+"]
        # myticks = myticks + [len(x) - 1]
        newlabels = newlabels[:-2]
        myticks = myticks[:-2]
        newlabels = newlabels + ["+"]
        myticks = myticks + [len(x) - 2]

    # if "problem" in list(x) and not "problem" in newlabels:
    #     newlabels = newlabels + ["bad"]
    #     myticks = myticks + [len(x) - 1]
    plt.xticks(myticks, newlabels, rotation = 85)
    plt.tick_params(labelsize=8)
    plt.title(sample)

    # print(_i, _j)
    # axs[_i][_j].set_xticks(myticks)
    # axs[_i][_j].set_xticklabels(newlabels)

    plt.savefig('../res/plot_count_cag/' + sample + '.png', dpi = 300)
    plt.close()
