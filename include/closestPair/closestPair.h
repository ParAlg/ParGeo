#pragma once

#include "parlay/sequence.h"
#include "pargeo/point.h"

namespace pargeo {

  template<int dim> std::pair<point<dim>, point<dim>>
  closestPair(parlay::sequence<point<dim>>& P, bool serial = false);

}
