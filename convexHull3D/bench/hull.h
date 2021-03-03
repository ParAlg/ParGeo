#ifndef HULL_H
#define HULL_H

#include "common/geometry.h"
#include "common/algebra.h"
#include "parlay/primitives.h"

typedef int intT;
typedef double floatT;

// CW-oriented 3d facet
struct facet3d {
  intT a, b, c;// Indices into P

  facet3d(intT aa, intT bb, intT cc, parlay::sequence<point<3>>& P): a(aa), b(bb), c(cc) {
    if (determinant3by3<floatT>(P[a], P[b], P[c]) > 0)
      swap(a, c);
  }

  facet3d(intT aa, intT bb, intT cc, parlay::slice<point<3>*, point<3>*>& P): a(aa), b(bb), c(cc) {
    if (determinant3by3<floatT>(P[a], P[b], P[c]) > 0)
      swap(a, c);
  }

  floatT area(parlay::sequence<point<3>>& P) {
    return crossProduct3d(P[b]-P[a], P[c]-P[a]).length()/2;
  }
};

parlay::sequence<facet3d> hull3d(parlay::sequence<point<3>> &S);

#endif
