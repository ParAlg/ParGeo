#pragma once

#include "geometry/point.h"
#include "geometry/algebra.h"
#include "common/algebra.h"
#include "parlay/primitives.h"

// CW-oriented 3d facet
template <class pt, class att>
struct _facet3d {
  size_t a, b, c;// Indices into P
  att attribute;

  _facet3d(size_t _a, size_t _b, size_t _c, parlay::slice<pt*, pt*> P): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(P[a], P[b], P[c]) > 0)
      swap(b, c);
  }
};

template <class pt>
using facet3d = _facet3d<pt, pargeo::_empty>;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P);
