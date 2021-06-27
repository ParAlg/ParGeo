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
        card = tokens[1]
        dataset = tokens1[1] + "_" + str(card)
        time = b['cpu_time']
        if dataset not in dataOut:
            dataOut[dataset] = dict()
        dataOut[dataset][method] = time

    return dataOut

def plotBar(data, title):
    print(data)
    # rows: methods
    # columns: data sets
    datasets = [k for k in data.keys()]
    methods = [m for m in data[datasets[0]]]
    bars = np.zeros((len(datasets), len(methods)),\
                    dtype = np.float64)
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            bars[i][j] = data[dataset][method]
    bars = bars.T

    X = np.arange(len(datasets))
    fig = plt.figure()

    w = 0.1
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red',\
              'tab:cyan', 'tab:purple', 'tab:brown', 'tab:pink',\
              'tab:gray', 'tab:olive']
    for i in range(len(methods)):
        print(bars[i])
        plt.bar(X + i * w, bars[i],\
                color = colors[i % len(methods)],\
                width = w)

    plt.xlabel("Data set")
    plt.ylabel("Time (ms)")

    plt.xticks(X, datasets)
    plt.legend(methods)
    plt.title(title)

    plt.show()

def main():
    file1 = ' '.join(sys.argv[1:])

    f1 = open(file1, "r")
    data1 = parse(json.load(f1))
    plotBar(data1, file1)

if __name__ == '__main__':
    main()
