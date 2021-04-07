#pragma once

#include "geometry/point.h"
#include "parlay/primitives.h" // ?
#include "parlay/sequence.h"

namespace sgInternal {
  struct _edge {
    size_t u,v;

    _edge(size_t _u, size_t _v): u(_u), v(_v) {
      if (u > v) std::swap(u,v);
    }

    _edge(): u(-1), v(-1) { }

    bool isEmpty() { return u == -1; }

    bool operator==(_edge e2) {
      return e2.u == u && e2.v == v;
    }

    bool operator!=(_edge e2) {
      return e2.u != u || e2.v != v;
    }
  };

}

using edge = sgInternal::_edge;

template<int dim>
parlay::sequence<edge> knnGraph(parlay::sequence<pargeo::point<dim>> &S, size_t k);

template<int dim>
parlay::sequence<edge> delaunayGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<edge> gabrielGraph(parlay::sequence<pargeo::point<dim>> &S);

template<int dim>
parlay::sequence<edge> betaSkeleton(parlay::sequence<pargeo::point<dim>> &S, double);
