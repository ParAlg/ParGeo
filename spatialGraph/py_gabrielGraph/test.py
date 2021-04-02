import numpy as np
from pygabriel import gabriel_graph
from datetime import datetime

start_time = None

def test(n, dim):
  A = np.random.random((n, dim))
  print("Data generated (python)")
  print("Gabriel graph, timer start")
  start_time = datetime.now()
  R = gabriel_graph(A)
  end_time = datetime.now()
  print('Duration: {}'.format(end_time - start_time))
  print("Result", R)

print("test 1 - 10")
test(10, 2)

print("\ntest 2 - 1m")
test(1000000, 2)
