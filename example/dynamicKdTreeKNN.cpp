#include <assert.h>
#include "dataset/uniform.h"
#include <iostream>
#include "pargeo/dynKdTree.h" /* Include dynamic kdTree */
#include "pargeo/getTime.h"
#include "parlay/parallel.h"
#include "pargeo/point.h"

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 4; // Data set dimensionality

  size_t n = 10000; // Number of data points

  size_t batches = 10; // Number of batch inserts

  /* Create a data set */

  auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);

  std::cout << "Test dynamic kdTree of dim " << dim << " on ";
  std::cout << "data of size " << n << "\n";

  /* Build initial dynamic kdTree by using half the points */

  pargeo::timer t;

  t.start();

  int base = n / 2;

  // unique_ptr<rootNode<dim, point<dim>>>
  //   tree(new rootNode<dim, point<dim>>(P, 0, base));

  pargeo::dynKdTree::rootNode<dim, pargeo::point<dim>> tree(P, 0, base);

  std::cout << "build-time = " << t.get_next() << "\n";

  /* Ration the remaining points into batches,
     and perform batch insertions */

  std::cout << "batches = " << batches << std::endl;

  int batchSize = (n - base) / batches;

  int inserted = 0;

  std::cout << "batch-size = " << batchSize << std::endl;

  for (int i = 0; i < batches; ++ i) {
    tree.insert(P, base + inserted, base + inserted + batchSize);
    inserted += batchSize;
  }

  std::cout << "size-after-insert = " << tree.size() << "\n";

  std::cout << "insert-time = " << t.get_next() << "\n";

  /* Perform 5-NN search for all n points */

  parlay::parallel_for (0, n, [&](size_t i) {
    tree.kNN(P[i], 5);
  });

  std::cout << "knn-time = " << t.get_next() << "\n";

  /* Erase the later half of the points from the tree */

  tree.erase(P, n / 2, n);

  std::cout << "size-after-erase = " << tree.size() << "\n";

  std::cout << "erase-time = " << t.get_next() << "\n";

  /* Compute 5-NN search for all n points again */

  parlay::parallel_for (0, n, [&](size_t i) {
    tree.kNN(P[i], 5);
  });

  std::cout << "knn-time = " << t.get_next() << "\n";

  assert(tree.check()); // Check if tree is still valid
  assert(tree.size() == n / 2);

}
