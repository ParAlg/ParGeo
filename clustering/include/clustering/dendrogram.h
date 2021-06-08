#pragma once

#include <tuple>
#include "parlay/sequence.h"
#include "pargeo/edge.h"
#include "pargeo/point.h"

namespace pargeo {

  typedef std::tuple<size_t, size_t, double, size_t> dendroNode;

  parlay::sequence<dendroNode> dendrogram(parlay::sequence<pargeo::wghEdge> &, size_t);

}
