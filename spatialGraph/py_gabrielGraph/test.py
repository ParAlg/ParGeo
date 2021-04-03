import numpy as np
from pygabriel import gabriel_graph
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

start_time = None

def test(n, name, savePlot = False):
  dim = 2
  print("\n" + name)

  points = np.random.random((n, dim))
  print("Data generated (python)")
  print("Gabriel graph, timer start")
  start_time = datetime.now()
  edges = gabriel_graph(points)
  end_time = datetime.now()
  print('Duration: {}'.format(end_time - start_time))
  print("#edges =", len(edges)/2)
  edges.shape = (-1, 2) # reshape to 2 x #edges

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

test(10, "test 1 - 10", False)
test(100, "test 2 - 100", True)
test(1000000, "test 3 - 100", False)
