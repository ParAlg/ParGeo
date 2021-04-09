import numpy as np
from pyknn import kdtknn
from datetime import datetime

start_time = None

def test(n, dim, k):
  A = np.random.random((n, dim))
  print("Data generated (python)")
  print("Knn, timer start")
  start_time = datetime.now()
  R = kdtknn(A, k)
  end_time = datetime.now()
  print('Duration: {}'.format(end_time - start_time))
  print("Result", R)

print("test 1 - 10, 2d")
test(10, 2, 2)

print("\ntest 2 - 10k, 2d")
test(10000, 3, 2)

print("\ntest 2 - 1m, 2d")
test(1000000, 3, 2)
