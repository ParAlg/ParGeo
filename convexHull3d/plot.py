from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import os.path

fig = plt.figure()
ax = plt.axes(projection='3d')

files = ["point.txt", "hull.txt", "other.txt"]
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
  if colors[I] == 'black':
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I], s=5);
  else:
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I]);

  if False:# len(labels) > 0:
    for p,label in enumerate(labels):
      ax.text(point1[p,0], point1[p,1], point1[p,2], label, size=7, zorder=1, color=colors[I])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("pargeo: my 3d hull")
plt.show()
