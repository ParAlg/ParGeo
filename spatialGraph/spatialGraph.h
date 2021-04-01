#pragma once

#include "geometry/point.h"
#include "parlay/primitives.h"

namespace sgInternal {
  struct _edge {
    size_t u,v;

    _edge(size_t _u, size_t _v): u(_u), v(_v) {
      if (u > v) swap(u,v);
    }

    _edge(): u(-1), v(-1) { }

    bool isEmpty() { return u < 0; }

    bool operator==(_edge e2) {
      return e2.u == u && e2.v == v;
    }

    bool operator!=(_edge e2) {
      return e2.u != u || e2.v != v;
    }
  };

  // std::ostream& operator<<(std::ostream& os, const _edge& e) {
  //   return os << "(" << e.u << "," << e.v << ")";
  // }

}

using edge = sgInternal::_edge;

// outputs edge list
template<int dim>
parlay::sequence<edge> spatialGraph(parlay::sequence<pargeo::point<dim>> &S);
