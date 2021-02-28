#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "knn.h"
#include "kdtKnn.h"

using namespace std;

template<int dim>
void knn(parlay::sequence<point<dim>> &P) {

  parlay::sequence<size_t> nnIdx = kdtKnn::kdtKnn<dim, point<dim>>(P, 3);
  parlay::sequence<size_t> nnIdx2 = kdtKnn::bruteforceKnn<dim, point<dim>>(P, 3);
  for(size_t i=0; i<10; ++i) {
    cout << P[i] << endl;
    for(size_t j=0; j<3; ++j)
      cout << nnIdx[i*3+j] << " ";
    cout << endl;
    for(size_t j=0; j<3; ++j)
      cout << nnIdx2[i*3+j] << " ";
    cout << endl;
  }
  cout << "..." << endl;

}

template void knn<2>(parlay::sequence<point<2>> &);
template void knn<3>(parlay::sequence<point<3>> &);
