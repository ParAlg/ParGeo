#include <assert.h>
#include "dataset/uniform.h"
#include <iostream>
#include "pargeo/kdTreeRange.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 4; // Data set dimensionality

  size_t n = 10000; // Number of data points

  /* Create a data set */

  auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);

  /* Build a tree */

  pargeo::kdNode<dim, pargeo::point<dim>>* tree =
    pargeo::buildKdTree<dim, pargeo::point<dim>>(P, true, false);

  // /* Spherical range query example
  //    surrounding P[0] with radius 0.1 */

  // parlay::sequence<size_t> elems1 = pargeo::kdTreeRange(P, tree, P[0], 0.1);

  // /* Rectangular range query example
  //    surrounding P[0] with half-length 0.1 */

  // parlay::sequence<size_t> elems2 = pargeo::kdTreeOrthRange(P, tree, P[0], 0.1);

  free(tree);
}
