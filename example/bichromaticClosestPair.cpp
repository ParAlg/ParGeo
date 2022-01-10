#include "dataset/uniform.h"
#include <iostream>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/kdTree.h" // Include kdTree
#include "pargeo/bccp.h" // Include BCCP
#include <tuple>

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 4; // Data set dimensionality

  size_t n = 10000; // Number of data points

  /* Create two point sets */

  auto red = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);
  auto blue = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);

  /* Build two kdTrees between the two point sets */

  pargeo::kdNode<dim, pargeo::point<dim>>* treeRed =
    pargeo::buildKdt<dim, pargeo::point<dim>>(red, true, false);

  pargeo::kdNode<dim, pargeo::point<dim>>* treeBlue =
    pargeo::buildKdt<dim, pargeo::point<dim>>(blue, true, true);

  /* Compute bichromatic closest pair between blue and red point sets */

  std::tuple<pargeo::point<dim>, pargeo::point<dim>, double> myBccp =
    bccp(treeRed, treeBlue);

  free(treeRed);
  free(treeBlue);
}
