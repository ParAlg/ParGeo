import matplotlib.pyplot as plt
import numpy as np
import json
import sys

def parse(data):
    bench = data['benchmarks']
    dataOut = dict()

    for b in bench:
        name = b['name']
        tokens = name.split("/")
        tokens1 = tokens[0].split("_")
        method = tokens1[0]
        dataset = "_".join(tokens1[1:])
        time = b['cpu_time']
        if dataset not in dataOut:
            dataOut[dataset] = dict()
        dataOut[dataset][method] = time

    return dataOut

def plotBar(data, title):
    datasets = [k for k in data.keys()]
    methods = [m for m in data[datasets[0]]]
    bars = np.zeros((len(datasets), len(methods)),\
                    dtype = np.float64)
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            bars[i][j] = data[dataset][method]

    ncol = int((len(datasets) / 2))
    fig, axes = plt.subplots(nrows=2, ncols=ncol, figsize=[15,7])

    axx = None
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',\
              'tab:cyan', 'tab:purple', 'tab:brown', 'tab:pink',\
              'tab:gray', 'tab:olive']
    for i, ax in enumerate(axes.flatten()):
        B = ax.bar(methods, bars[i])
        for j, bar in enumerate(B):
            bar.set_label(methods[j])
            bar.set_color(colors[j])
        ax.set(title=datasets[i], ylabel='time (ms)')
        ax.tick_params(labelbottom = False, bottom = False)
        axx = ax
    handles, labels = axx.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=len(methods),\
               fontsize=12, bbox_to_anchor=(0.5, 0.95))

    fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top = 0.85)
    #plt.show()
    plt.savefig(title + ".png")

def main():
    file1 = ' '.join(sys.argv[1:])

    f1 = open(file1, "r")
    data1 = parse(json.load(f1))
    plotBar(data1, file1)

if __name__ == '__main__':
    main()
