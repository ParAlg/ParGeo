#pragma once

#include "geometry/point.h"
#include "geometry/algebra.h"
#include "parlay/primitives.h"

// CW-oriented 3d facet
template <class pt, class att>
struct _facet3d {
  using pointT = pt;
  pt a, b, c;// Indices into P
  att attribute;

  _facet3d(pt _a, pt _b, pt _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > 0)
      std::swap(b, c);
  }
};

template <class pt>
using facet3d = _facet3d<pt, pargeo::_empty>;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P);

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P, size_t);
