#pragma once

#include "point.h"

namespace pargeo {
  // All inputs are row-major

  template<class pt>
  inline typename pt::floatT determinant2by2(pt a, pt b) {
    return a[0]*b[1] - a[1]*b[0];
  }

  template<class pt>
  inline void inverse2by2(pt a, pt b) {
    typename pt::floatT c = 1/determinant2by2(a, b);
    typename pt::floatT m[4];
    m[0] =  b[1]*c; m[1] = -a[1]*c;
    m[2] = -b[0]*c; m[3] =  a[0]*c;
    a[0] = m[0]; a[1] = m[1];
    b[0] = m[2]; b[1] = m[3];
  }

  template <class pt>
  inline pt crossProduct3d(pt p, pt q) {
    pt r;
    r[0] = p[1]*q[2] - p[2]*q[1];
    r[1] = p[2]*q[0] - p[0]*q[2];
    r[2] = p[0]*q[1] - p[1]*q[0];
    return r;
  }

  template<class pt>
  inline typename pt::floatT determinant3by3(pt a, pt b, pt c) {
    return a[0] * (b[1]*c[2] - c[1]*b[2]) -
      a[1] * (b[0]*c[2] - b[2]*c[0]) +
      a[2] * (b[0]*c[1] - b[1]*c[0]);
  }

  template<class pt>
  inline typename pt::floatT determinant4by4(pt a, pt b, pt c, pt d) {
    return
      a[3] * b[2] * c[1] * d[0] - a[2] * b[3] * c[1] * d[0] -
      a[3] * b[1] * c[2] * d[0] + a[1] * b[3] * c[2] * d[0] +
      a[2] * b[1] * c[3] * d[0] - a[1] * b[2] * c[3] * d[0] -
      a[3] * b[2] * c[0] * d[1] + a[2] * b[3] * c[0] * d[1] +
      a[3] * b[0] * c[2] * d[1] - a[0] * b[3] * c[2] * d[1] -
      a[2] * b[0] * c[3] * d[1] + a[0] * b[2] * c[3] * d[1] +
      a[3] * b[1] * c[0] * d[2] - a[1] * b[3] * c[0] * d[2] -
      a[3] * b[0] * c[1] * d[2] + a[0] * b[3] * c[1] * d[2] +
      a[1] * b[0] * c[3] * d[2] - a[0] * b[1] * c[3] * d[2] -
      a[2] * b[1] * c[0] * d[3] + a[1] * b[2] * c[0] * d[3] +
      a[2] * b[0] * c[1] * d[3] - a[0] * b[2] * c[1] * d[3] -
      a[1] * b[0] * c[2] * d[3] + a[0] * b[1] * c[2] * d[3];
  }

} // End namespace
