#pragma once

#include "geometry/point.h"
#include "geometry/algebra.h"
#include "common/algebra.h"
#include "parlay/primitives.h"

// CW-oriented 3d facet
template <class pt>
struct facet3d {
  size_t a, b, c;// Indices into P

  facet3d(size_t _a, size_t _b, size_t _c, parlay::slice<pt*, pt*> P): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(P[a], P[b], P[c]) > 0)
      swap(a, c);
  }
};

parlay::sequence<facet3d<pargeo::point<3>>> hull3d(parlay::sequence<pargeo::point<3>> &P);
