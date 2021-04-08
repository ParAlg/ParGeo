#include "geometry/point.h"
#include "parlay/primitives.h"

namespace pargeo {

  template<int dim>
  parlay::sequence<size_t> knn(parlay::sequence<pargeo::point<dim>> &S, size_t k);

}
