import numpy as np
from pypargeo import loadPoints
from pypargeo import HDBSCAN
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from sklearn.datasets import make_blobs

points, y = make_blobs(n_samples=20, centers=3, n_features=2, random_state=0)

dendro = HDBSCAN(points, 3)

fig = plt.figure(figsize=(3,2.3))

dn = dendrogram(dendro, color_threshold=0, no_labels=True)

plt.show()

fig.savefig("hdbscan.pdf")
