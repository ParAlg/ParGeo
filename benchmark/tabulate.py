import matplotlib.pyplot as plt
import numpy as np
import json
import sys

def parse(data):
    bench = data['benchmarks']
    dataOut = dict()

    for b in bench:
        method = b['name']
        time = b['real_time']
        dataOut[method] = time

    return dataOut

def main():
    files = ' '.join(sys.argv[1:])
    tokens = files.split(" ")

    if len(tokens) <= 1:
        print("need two files, \"1-thread.json multi-thread.json\"")
        exit(1)

    file1 = tokens[0]
    file2 = tokens[1]

    f1 = open(file1, "r")
    data1 = parse(json.load(f1))
    f2 = open(file2, "r")
    data2 = parse(json.load(f2))

    print(files)

    methods = list(data1.keys())

    print("\nmethods ---------")
    print(", ".join(methods))
    print("")

    for method in methods:
        su = data1[method] / data2[method]
        print(method + " T1 {:.2f} Tm {:.2f} Su {:.2f}x".format(data1[method], data2[method], su))

if __name__ == '__main__':
    main()
