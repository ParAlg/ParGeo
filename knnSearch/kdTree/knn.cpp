#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "geometry/point.h"
#include "knn.h"
#include "kdtKnn.h"

using namespace std;

template<int dim>
parlay::sequence<size_t> knn(parlay::sequence<pargeo::point<dim>> &P, size_t k) {
  if (k < 2 || k > P.size()) {
    cout << "Error, k = " << k << " k must range from 2 to #-data-points, abort" << endl;
    abort();
  }

  parlay::sequence<size_t> nnIdx = kdtKnn::kdtKnn<dim, pargeo::point<dim>>(P, k);
  //parlay::sequence<size_t> nnIdx2 = kdtKnn::bruteforceKnn<dim, pargeo::point<dim>>(P, k);

  return nnIdx;
}

template parlay::sequence<size_t> knn<2>(parlay::sequence<pargeo::point<2>> &, size_t);
template parlay::sequence<size_t> knn<3>(parlay::sequence<pargeo::point<3>> &, size_t);
