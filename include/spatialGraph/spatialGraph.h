#pragma once

#include "parlay/sequence.h"
#include "pargeo/edge.h"
#include "pargeo/point.h"

namespace pargeo {

template<int dim>
parlay::sequence<pargeo::dirEdge> knnGraph(parlay::sequence<pargeo::point<dim>> &S, size_t k);

template<int dim>
parlay::sequence<pargeo::edge> delaunayGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<pargeo::edge> gabrielGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<pargeo::edge> betaSkeleton(parlay::sequence<pargeo::point<dim>> &S, double);

template<int dim>
parlay::sequence<pargeo::edge> spanner(parlay::sequence<pargeo::point<dim>> &S, double);

}
