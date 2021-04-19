import numpy as np
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from geograph import delaunay_graph
from geograph import gabriel_graph
from geograph import knn_graph
from geograph import beta_skeleton

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
test(points, "rand-delaunay-100", delaunay_graph, True)
test(points, "rand-gabriel-100", gabriel_graph, True)
test(points, "rand-knn1-100", knn_graph, True, 1)
test(points, "rand-knn2-100", knn_graph, True, 2)
test(points, "rand-knn10-100", knn_graph, True, 10)
test(points, "beta-skeleton0.5-100", beta_skeleton, True, 0.5)
test(points, "beta-skeleton1-100", beta_skeleton, True, 1)
test(points, "beta-skeleton2-100", beta_skeleton, True, 2)
test(points, "beta-skeleton5-100", beta_skeleton, True, 5)
