def plot_repeat_dist(tsv, output_dir, sample_name, N):
    """
    tsv:
        col1 col2
        repeat_length   count
    N = 5  # 1 tick every 5
    """
    # output plot:
    output_plot = os.path.join(output_dir, sample_name + ".png")
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)

    # We can set the number of bins with the *bins* keyword argument.
    df = pd.read_csv(tsv, sep = "\t")
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
    n, bins, patches = plt.hist(x, weights = weight, bins=len(x), range = (0, len(x)))

    # sparse the x-axis ticks:
    myticks = [i for i in range(len(x)) if not i%N]
    newlabels = [i for i in range(len(x)) if not i%N]

    if "plus" in list(x) and not "plus" in newlabels:
        newlabels = newlabels[:-2]
        myticks = myticks[:-2]
        newlabels = newlabels + ["+"]
        myticks = myticks + [len(x) - 2]

    plt.xticks(myticks, newlabels, rotation = 85)
    plt.tick_params(labelsize=8)
    plt.title(sample_name)

    plt.savefig(output_plot, dpi = 300)
    plt.close()

def print_test():
    return("inside here")
