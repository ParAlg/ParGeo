#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "geometry/point.h"
#include "knn.h"
#include "kdtKnn.h"

using namespace std;

template<int dim>
parlay::sequence<size_t> pargeo::knn(parlay::sequence<pargeo::point<dim>> &P, size_t k) {
  if (k < 1 || k > P.size()-1) {
    cout << "Error, k = " << k << " k must range from 1 to #-data-points-1, abort" << endl;
    abort();
  }

  parlay::sequence<size_t> nnIdx = kdtKnn::kdtKnn<dim, pargeo::point<dim>>(P, k+1);
  //parlay::sequence<size_t> nnIdx2 = kdtKnn::bruteforceKnn<dim, pargeo::point<dim>>(P, k+1);

  return nnIdx;
}

template parlay::sequence<size_t> pargeo::knn<2>(parlay::sequence<pargeo::point<2>> &, size_t);
template parlay::sequence<size_t> pargeo::knn<3>(parlay::sequence<pargeo::point<3>> &, size_t);
