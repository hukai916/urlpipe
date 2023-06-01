import matplotlib.pyplot as plt
import os
import pandas as pd
from Bio.Seq import Seq

def plot_repeat_dist(csv, output_file, sample_name, N, bin_number = 250):
    """
    csv:
        col1 col2
        repeat_length   count
    N = 5  # 1 tick every 5
    """
    # output plot:
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # read in stat csv: ensure 'plus' and 'problem' are there
    df = pd.read_csv(csv, sep = ",", header = None)
    x, y = [], []
    _dict = {"plus": 0, "problem": 0} # ensure that 'plus' and 'problem' are there
    print(df.iloc[:,0])
    for i in range(len(df.iloc[:,0])):
        print(df.iloc[i, 0])
        if not df.iloc[i, 0] in ["plus", "problem"]:
            _dict[str(round(float((df.iloc[i, 0]))))] = df.iloc[i, 1] # round the float to int, for plotting purpose below
        else:
            _dict[str(df.iloc[i, 0])] = df.iloc[i, 1]

    # determine max count:
    max_count = max([int(i) for i in df.iloc[:,0] if str(i).replace(".", "").isdigit()])
    # plus, problem = 0, 0
    # for i, v in enumerate(df.iloc[:, 0]):
    #     if v == "plus":
    #         plus = df.iloc[i, 1]
    #     elif v == "problem":
    #         problem = df.iloc[i, 1]

    # determine bin_n considering 'auto' option:
    if bin_number == 'auto': # if 'auto', use max_count as bin_number
        bin_n = max_count # include "plus" and "problem"
    else:
        bin_n = int(bin_number)
        assert max_count <= bin_n, "X-axis scale too small!"

    # determine x and weight (y), add N number of gap between count and "plus" and "problem" for better visualization
    x, y = [], []
    for i, v in enumerate(list(range(bin_n + N)) + ["plus"] + [0] * N +  ["problem"]):
        x = x + [str(i)]
        if not v in ["plus", "problem"]:
            if not str(i) in _dict:
                y = y + [0]
            else:
                y = y + [_dict[str(i)]]
        elif v == "plus":
            y = y + [_dict["plus"]]
        elif v == "problem":
            y = y + [_dict["problem"]]

    # sparse the ticks:
    myticks = [i for i in range(bin_n) if not i%N]
    mylabels = [i for i in range(bin_n) if not i%N]
    if not "plus" in mylabels:
        mylabels = mylabels + ["plus"]
        myticks = myticks + [bin_n + N]
    if not "problem" in mylabels:
        mylabels = mylabels + ["problem"]
        myticks = myticks + [bin_n + 2* N]

    # plot:
    n, bins, patches = plt.hist(x, weights = y, bins = bin_n + 2 * N, range = (0, bin_n + 2 * N))
    plt.xticks(myticks, mylabels, rotation = 85)
    plt.tick_params(labelsize = 8)

    plt.savefig(output_file, dpi = 300, bbox_inches = "tight")
    plt.close()

def print_test():
    return("inside here")

def reverse_complement(seq):
    dna = Seq(seq)
    return str(dna.reverse_complement())

def get_sample_name(x, prefix, suffix):
    if prefix == "" and suffix == "":
        return(x)
    elif prefix == "":
        return(x.split(suffix)[0])
    elif suffix == "":
        return(x.split(prefix)[1])
    else:
        return(x.split(prefix)[1].split(suffix)[0])

# def get_std(x):
#     """
#     Return weighted std. If empty input, return na.
#     """
#     try:
#         _average = np.average((x - np.average(x, weights = df.loc[x.index, sample_name])) ** 2, weights = df.loc[x.index, sample_name])
#     except:
#         _average = np.nan
#     return(np.sqrt(_average))

def indel_filter(indel, no_indel, indel_cutoff = 0.5, add = False):
    """
    For indel reads belong to UMI groups that also have reads in no_indel:
    if indel reads percentage is above 50% (cutoff > 0.5), move them out from the indel category.
    Both indel and no_indel are lists containing Bio.Seq record objects.

    If add = True, will add moved records from indel to no_indel and return [indel_res, no_indel_res]; otherwise return indel_res only.
    """
    
    # keep track of reads per umi:
    dict_umi = {} # [reads_per_umi_no_indel, reads_per_umi_indel]
    for record in no_indel:
        umi = record.name.split("_")[1]
        if not umi in dict_umi:
            dict_umi[umi] = [1, 0]
        else:
            dict_umi[umi][0] += 1
    for record in indel:
        umi = record.name.split("_")[1]
        if not umi in dict_umi:
            dict_umi[umi] = [0, 1]
        else:
            dict_umi[umi][1] += 1
    
    indel_res = []
    for record in indel:
        umi = record.name.split("_")[1]
        if dict_umi[umi][1] / (dict_umi[umi][0] + dict_umi[umi][1]) > float(indel_cutoff):
            indel_res.append(record)
        if add:
            no_indel.append(record)
    if not add:
        return(indel_res)
    else:
        return([indel_res, no_indel])
