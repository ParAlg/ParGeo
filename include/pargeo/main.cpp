#include <assert.h>
#include <iostream>
#include <memory>
#include <vector>

#include "dynKdTree.h"
#include "getTime.h"

template<int dim>
class point: public dynKdTree::coordinate<dim> {

public:

  point(double* _data): dynKdTree::coordinate<dim>(_data) { }

};


int main() {

  using namespace std;

  using namespace dynKdTree;

  static const int dim = 2;

  int n = 2000000;

  /* Generate test data */

  unique_ptr<double[]> rawData(new double[n]);

  for (int i = 0; i < n; ++ i) {
    rawData[i] = rand() / (double)RAND_MAX;
  }

  container<point<dim>> points;

  for (int i = 0; i < n / dim; ++ i) {
    points.push_back(point<dim>(&rawData[dim * i]));
  }

  /* Insert to the kdtree */

  pargeo::timer t;
  t.start();

  unique_ptr<splitNode<dim, point<dim>>>
    tree(new splitNode<dim, point<dim>>(points, 0, n / 4));

  std::cout << "build-time = " << t.get_next() << "\n";

  tree->insert(points, n / 4, n / 2);

  std::cout << "size-after-insert = " << tree->size() << "\n";

  std::cout << "insert-time = " << t.get_next() << "\n";

  tree->erase(points, n / 4, n / 2);

  std::cout << "size-after-erase = " << tree->size() << "\n";

  std::cout << "erase-time = " << t.get_next() << "\n";

  assert(tree->size() == n / 4);

  return 0;

}
