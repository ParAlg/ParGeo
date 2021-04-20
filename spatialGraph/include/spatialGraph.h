#pragma once

#include "geometry/point.h"
#include "parlay/primitives.h" // ?
#include "parlay/sequence.h"
#include "geometry/edge.h"

namespace pargeo {

using edge = pargeo::_edge<false, unsigned int>;
using dirEdge = pargeo::_edge<true, unsigned int>;

template<int dim>
parlay::sequence<pargeo::dirEdge> knnGraph(parlay::sequence<pargeo::point<dim>> &S, size_t k);

template<int dim>
parlay::sequence<pargeo::edge> delaunayGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<pargeo::edge> gabrielGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<pargeo::edge> betaSkeleton(parlay::sequence<pargeo::point<dim>> &S, double);

}
