#include "common/geometry.h"
#include "parlay/primitives.h"

template<int dim>
parlay::sequence<size_t> knn(parlay::sequence<point<dim>> &S, size_t k);
