// Pargeo
// todo deprecated

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "geometry.h"

inline point<3> crossProduct3d(point<3> p, point<3> q) {
  point<3> r;
  r[0] = p[1]*q[2] - p[2]*q[1];
  r[1] = p[2]*q[0] - p[0]*q[2];
  r[2] = p[0]*q[1] - p[1]*q[0];
  return r;
}

//takes row-major 2by2 matrix m
template<class floatT>
inline floatT determinant2by2(floatT* m) {
  return m[0]*m[3] - m[1]*m[2];
}

//takes row-major 2by2 matrix m, output replaces m
template<class floatT>
inline void inverse2by2(floatT* m) {
  auto coeff = 1/determinant2by2(m);
  floatT mInv[4];
  mInv[0] = m[3]*coeff;
  mInv[1] = -m[1]*coeff;
  mInv[2] = -m[2]*coeff;
  mInv[3] = m[0]*coeff;
  m[0] = mInv[0]; m[1] = mInv[1]; m[2] = mInv[2]; m[3] = mInv[3];
}

//takes row-major 3by3 matrix m
template<class floatT>
inline floatT determinant3by3(floatT* m) {
  return m[0*3+0] * (m[1*3+1] * m[2*3+2] - m[2*3+1] * m[1*3+2]) -
    m[0*3+1] * (m[1*3+0] * m[2*3+2] - m[1*3+2] * m[2*3+0]) +
    m[0*3+2] * (m[1*3+0] * m[2*3+1] - m[1*3+1] * m[2*3+0]);
}

//takes rows of matrix m as points
template<class floatT>
inline floatT determinant3by3(point<3> a, point<3> b, point<3> c) {
  return a[0] * (b[1] * c[2] - c[1] * b[2]) -
    a[1] * (b[0] * c[2] - b[2] * c[0]) +
    a[2] * (b[0] * c[1] - b[1] * c[0]);
}

//takes row-major 3by3 matrix m, output replaces m
template<class floatT>
inline void inverse3by3(floatT* m) {
  floatT invDet = 1/determinant3by3(m);
  floatT mInv[9];
  mInv[0*3+0] = (m[1*3+1] * m[2*3+2] - m[2*3+1] * m[1*3+2]) * invDet;
  mInv[0*3+1] = (m[0*3+2] * m[2*3+1] - m[0*3+1] * m[2*3+2]) * invDet;
  mInv[0*3+2] = (m[0*3+1] * m[1*3+2] - m[0*3+2] * m[1*3+1]) * invDet;
  mInv[1*3+0] = (m[1*3+2] * m[2*3+0] - m[1*3+0] * m[2*3+2]) * invDet;
  mInv[1*3+1] = (m[0*3+0] * m[2*3+2] - m[0*3+2] * m[2*3+0]) * invDet;
  mInv[1*3+2] = (m[1*3+0] * m[0*3+2] - m[0*3+0] * m[1*3+2]) * invDet;
  mInv[2*3+0] = (m[1*3+0] * m[2*3+1] - m[2*3+0] * m[1*3+1]) * invDet;
  mInv[2*3+1] = (m[2*3+0] * m[0*3+1] - m[0*3+0] * m[2*3+1]) * invDet;
  mInv[2*3+2] = (m[0*3+0] * m[1*3+1] - m[1*3+0] * m[0*3+1]) * invDet;
  for(int i=0; i<9; ++i) m[i] = mInv[i];
}

//takes row-major 4by4 matrix m
template<class floatT>
inline floatT determinant4by4(floatT* m) {
  return
    m[0*4+3] * m[1*4+2] * m[2*4+1] * m[3*4+0] - m[0*4+2] * m[1*4+3] * m[2*4+1] * m[3*4+0] -
    m[0*4+3] * m[1*4+1] * m[2*4+2] * m[3*4+0] + m[0*4+1] * m[1*4+3] * m[2*4+2] * m[3*4+0] +
    m[0*4+2] * m[1*4+1] * m[2*4+3] * m[3*4+0] - m[0*4+1] * m[1*4+2] * m[2*4+3] * m[3*4+0] -
    m[0*4+3] * m[1*4+2] * m[2*4+0] * m[3*4+1] + m[0*4+2] * m[1*4+3] * m[2*4+0] * m[3*4+1] +
    m[0*4+3] * m[1*4+0] * m[2*4+2] * m[3*4+1] - m[0*4+0] * m[1*4+3] * m[2*4+2] * m[3*4+1] -
    m[0*4+2] * m[1*4+0] * m[2*4+3] * m[3*4+1] + m[0*4+0] * m[1*4+2] * m[2*4+3] * m[3*4+1] +
    m[0*4+3] * m[1*4+1] * m[2*4+0] * m[3*4+2] - m[0*4+1] * m[1*4+3] * m[2*4+0] * m[3*4+2] -
    m[0*4+3] * m[1*4+0] * m[2*4+1] * m[3*4+2] + m[0*4+0] * m[1*4+3] * m[2*4+1] * m[3*4+2] +
    m[0*4+1] * m[1*4+0] * m[2*4+3] * m[3*4+2] - m[0*4+0] * m[1*4+1] * m[2*4+3] * m[3*4+2] -
    m[0*4+2] * m[1*4+1] * m[2*4+0] * m[3*4+3] + m[0*4+1] * m[1*4+2] * m[2*4+0] * m[3*4+3] +
    m[0*4+2] * m[1*4+0] * m[2*4+1] * m[3*4+3] - m[0*4+0] * m[1*4+2] * m[2*4+1] * m[3*4+3] -
    m[0*4+1] * m[1*4+0] * m[2*4+2] * m[3*4+3] + m[0*4+0] * m[1*4+1] * m[2*4+2] * m[3*4+3];
}

#endif
