import numpy as np
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from pypargeo import DelaunayGraph
from pypargeo import GabrielGraph
from pypargeo import KnnGraph
from pypargeo import BetaSkeleton

start_time = None

def test(points, name, gen, savePlot = False, k = -1):
  dim = 2
  print("\n" + name)

  start_time = datetime.now()
  if k < 0:
    edges = gen(points)
  else:
    edges = gen(points, k)
  end_time = datetime.now()
  print("graph-gen-time: {}".format(end_time - start_time))
  print("#-edges =", len(edges)/2)

  if savePlot:
    lc = LineCollection(points[edges])
    fig = plt.figure()
    plt.gca().add_collection(lc)
    plt.xlim(points[:,0].min()-.1, points[:,0].max()+.1)
    plt.ylim(points[:,1].min()-.1, points[:,1].max()+.1)
    plt.plot(points[:,0], points[:,1], 'ro', markersize=2)
    # labels = [str(i) for i in range(points.shape[0])]
    # for p,label in enumerate(labels):
    #   plt.text(points[p,0], points[p,1], label, size=12)
    plt.gca().set_aspect('equal', adjustable='box')
    fig.savefig(name+".png")

points = np.random.random((100, 2))
test(points, "rand-delaunay-100", DelaunayGraph, True)
test(points, "rand-gabriel-100", GabrielGraph, True)
test(points, "rand-knn1-100", KnnGraph, True, 1)
test(points, "rand-knn2-100", KnnGraph, True, 2)
test(points, "rand-knn10-100", KnnGraph, True, 10)
test(points, "beta-skeleton0.5-100", BetaSkeleton, True, 0.5)
test(points, "beta-skeleton1-100", BetaSkeleton, True, 1)
test(points, "beta-skeleton2-100", BetaSkeleton, True, 2)
test(points, "beta-skeleton5-100", BetaSkeleton, True, 5)
