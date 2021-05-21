#pragma once

#include "parlay/sequence.h"
#include "pargeo/edge.h"
#include "pargeo/point.h"

namespace pargeo {

  // todo introduce weight into edge
  template<int dim>
  parlay::sequence<pargeo::edge> hdbscan(parlay::sequence<pargeo::point<dim>> &, size_t);

}
