#pragma once

#include "pargeo/point.h"
#include "parlay/utilities.h"
#include "parlay/sequence.h"

namespace pargeo {

  template<int dim>
  parlay::sequence<pargeo::point<dim>> uniformInPolyPoints(size_t n, size_t shape, double scale = 1.0);

  template<int dim>
  parlay::sequence<pargeo::point<dim>> uniformOnPolyPoints(size_t n, size_t shape, double thickness, double scale = 1.0);

}
