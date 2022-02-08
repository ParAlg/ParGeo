from mpl_toolkits import mplot3d
import matplotlib.colors as colors
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import os.path

fig = plt.figure(figsize=(30,20))
ax = plt.axes(projection='3d')

files = ["point.txt", "red.txt", "blue.txt"]
colors = ['black', 'red', 'blue']

for I,file in enumerate(files):
  if not os.path.isfile(file):
    continue
  pointFile = open(file)
  labels = list()
  point1 = list()
  for i,line in enumerate(pointFile):
      tokens = line.split(' ')
      if len(tokens) < 3:
          continue
      point1.append([float(tokens[0]), float(tokens[1]), float(tokens[2])]);
      if len(tokens) >= 4:
        labels.append(tokens[3][:-1])
  pointFile.close()

  point1 = np.asarray(point1)
  if point1.shape[0] <= 0:
      continue
  if colors[I] == 'blue':
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I], s=20, alpha=0.3);
  if colors[I] == 'black':
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I], s=3);
  else:
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I], s=50, alpha=0.3);

  if len(labels) > 0:
    for p,label in enumerate(labels):
      ax.text(point1[p,0], point1[p,1], point1[p,2], label, size=8, zorder=1, color=colors[I])

pointFile = open("facet.txt")
labels = list()
vertices = list()
for i,line in enumerate(pointFile):
  tokens = line.split(' ')
  if len(tokens) < 3:
    continue
  vertices.append([float(tokens[0]), float(tokens[1]), float(tokens[2])]);
pointFile.close()

facets = list()
for i in range(int(len(vertices)/3)):
  vtx = np.zeros((3,3))
  vtx[0] = np.array(vertices[i*3])
  vtx[1] = np.array(vertices[i*3+1])
  vtx[2] = np.array(vertices[i*3+2])
  facets.append(vtx)
tri = mplot3d.art3d.Poly3DCollection(facets, alpha=0.1)
tri.set_color('yellow')
tri.set_edgecolor('k')
ax.add_collection3d(tri)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("pargeo: my 3d hull")
plt.show()
#plt.savefig("hull.png")
