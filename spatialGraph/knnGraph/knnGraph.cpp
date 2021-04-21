#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "pargeo/point.h"
#include "pargeo/kdtKnn.h"
#include "spatialGraph/spatialGraph.h"

using namespace std;

template<int dim>
parlay::sequence<pargeo::dirEdge> pargeo::knnGraph(parlay::sequence<pargeo::point<dim>> &P, size_t k) {
  using namespace parlay;

  if (k < 1 || k > P.size()-1)
    throw std::runtime_error("k must range from 1 to n-1");

  // Call knn
  parlay::sequence<size_t> nnIdx = kdtKnn::kdtKnn<dim, pargeo::point<dim>>(P, k+1);

  // Convert to edge list
  auto edges = parlay::sequence<dirEdge>(k * P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      size_t jj = 0;
			      for(size_t j=0; j<k+1; ++j) {
				if (nnIdx[i*(k+1)+j] != i) {
				  edges[i*k + jj] = dirEdge(i, nnIdx[i*(k+1)+j]);
				  jj++;
				}
			      }
			    });

  return edges;
}

template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<2>(parlay::sequence<pargeo::point<2>> &, size_t);
template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<3>(parlay::sequence<pargeo::point<3>> &, size_t);
template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<4>(parlay::sequence<pargeo::point<4>> &, size_t);
template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<5>(parlay::sequence<pargeo::point<5>> &, size_t);
template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<6>(parlay::sequence<pargeo::point<6>> &, size_t);
template parlay::sequence<pargeo::dirEdge> pargeo::knnGraph<7>(parlay::sequence<pargeo::point<7>> &, size_t);
