#pragma once

#include "parlay/sequence.h"
#include "pargeo/edge.h"
#include "pargeo/point.h"

namespace pargeo {

  template<int dim>
  parlay::sequence<pargeo::wghEdge> hdbscan(parlay::sequence<pargeo::point<dim>> &, size_t);

}
