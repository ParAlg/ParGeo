import numpy as np
import matplotlib.pyplot as plt

class point:
  def __init__(self, x, y, label="", color="k"):
    self.x = x
    self.y = y
    self.label = label
    self.color = color
  def __str__(self):
    return "(" + str(self.x) + ", " + str(self.y) + ")"

class circle:
  def __init__(self, center, radius, label="", color="k"):
    self.center = center
    self.radius = radius
    self.label = label
    self.color = color
  def __str__(self):
    return "(" + str(self.center) + ", rad=" + str(self.radius) + ")"

class snapshot:
    path = "./plots/"
    xMin = -3
    xMax = 3
    yMin = -3
    yMax = 3

    def __init__(self, points, circle, idNum):
        self.points = points;
        self.circle = circle;
        self.idNum = idNum;

    def plot(self):
      fig, ax = plt.subplots(figsize=(10,10))

      for pt in self.points:
        plt.plot([pt.x], [pt.y], marker='o', markersize=3, color=pt.color)
        plt.text(pt.x*1.02, pt.y*1.02, str(pt.label), fontsize=10)

      ax.add_patch(plt.Circle(np.array([self.circle.center.x, self.circle.center.y]), radius=self.circle.radius, edgecolor='k', facecolor='w', alpha=0.5))

      #scale_factor = 3
      #xmin, xmax = plt.xlim()
      #ymin, ymax = plt.ylim()
      #plt.xlim(xmin * scale_factor - xmin, xmax * scale_factor - xmax)
      #plt.ylim(ymin * scale_factor - ymin, ymax * scale_factor - ymax)
      plt.xlim(self.xMin, self.xMax)
      plt.ylim(self.yMin, self.yMax)

      plt.grid()
      #plt.show()
      plt.savefig(self.path + str(self.idNum) + ".png", format="png")
      plt.close()

def enumerateSnap():
  path = "./"
  snaps = list()
  import os
  for file in os.listdir(path):
    if file.endswith(".txt"):
      snaps.append(os.path.join(path, file))
  #print(snaps)
  snapss = list()
  for i in range(len(snaps)):#hack
    snapss.append(os.path.join(path, str(i)+".txt"))
  return snapss

def parseSnap(path, idNum):
  f = open(path)
  lines = f.readlines()

  #circle
  tokens = lines[0].split();
  c = circle(
      point(float(tokens[0]), float(tokens[1])),
      float(tokens[2]))

  tokens = lines[1].split();
  ci = int(tokens[0])

  #points
  points = list()
  for line in lines[2:]:
    line = line[:-1]
    tokens = line.split(" ")
    idx = int(tokens[2])
    points.append(point(
        float(tokens[0]),
        float(tokens[1]),
        str(idx),
        "r" if idx==ci else "k"
      ))
  return snapshot(points, c, idNum)

files = enumerateSnap()

for i,f in enumerate(files):
  print("snapshot " + str(i))
  parseSnap(f,i).plot()
