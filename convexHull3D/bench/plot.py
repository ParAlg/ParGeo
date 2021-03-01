from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection='3d')

files = ["point.txt", "hull.txt"]
colors = ['blue', 'red']

for I,file in enumerate(files):
  pointFile = open(file)
  point1 = list()
  for i,line in enumerate(pointFile):
      tokens = line.split(' ')
      if len(tokens) < 3:
          continue
      point1.append([float(tokens[0]), float(tokens[1]), float(tokens[2])]);
  pointFile.close()
  point1 = np.asarray(point1)
  if point1.shape[0] <= 0:
      continue
  if colors[I] == 'black':
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I], s=5);
  else:
    ax.scatter3D(point1[:,0].T, point1[:,1].T, point1[:,2].T, color=colors[I]);

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.title("pargeo: my 3d hull")
plt.show()
