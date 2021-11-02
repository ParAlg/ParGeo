#include <iostream>
#include <memory>
#include <vector>

#include "dynKdTree.h"

template<int dim>
class point: public dynKdTree::coordinate<dim> {

public:

  point(double* _data): dynKdTree::coordinate<dim>(_data) { }

};


int main() {

  using namespace std;

  using namespace dynKdTree;

  static const int dim = 2;

  int n = 80;

  /* Generate test data */

  unique_ptr<double[]> rawData(new double[n]);

  for (int i = 0; i < n; ++ i) {
    rawData[i] = rand() / (double)RAND_MAX;
  }

  vector<point<dim>> points;

  for (int i = 0; i < n / dim; ++ i) {
    points.push_back(point<dim>(&rawData[dim * i]));
  }

  /* Insert to the kdtree */

  // boundingBox<dim> bb(points);
  // cout << bb.topLeft[0] << "," << bb.lowerRight[0] << endl;
  // cout << bb.topLeft[1] << "," << bb.lowerRight[1] << endl;

  unique_ptr<splitNode<dim, point<dim>>>
    tree(new splitNode<dim, point<dim>>(points, 0, n / 4));

  tree->insert(points, n / 4, n / 2);

  return 0;
}
