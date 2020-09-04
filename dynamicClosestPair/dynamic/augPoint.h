// This code is part of the project "A Parallel Batch-Dynamic Data Structure
// for the Closest Pair Problem"
// Copyright (c) 2020 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef AUG_POINT_H
#define AUG_POINT_H

#include "grid.h"
#include "geometry.h"

/**
 *  An augmented point class.
 */
template<int dim>
struct augPoint {
  typedef double floatT;
  typedef point<dim> geoPointT;
  typedef augPoint<dim> pointT;
  static const intT clearFlag = 0;
  static const intT emptyFlag = -1;
  geoPointT p;

  augPoint(geoPointT pp): p(pp) {};
  augPoint() {setEmpty();};

  inline bool isEmpty() {return p.isEmpty();}
  inline void setEmpty() {p = geoPointT();}

  floatT operator[](intT i) {return p[i];}
  floatT pointDist(pointT q) {return p.pointDist(q.p);}
  floatT pointDist(geoPointT q) {return p.pointDist(q);}
  geoPointT pt() {return p;}
  floatT* x() {return p.x;}
  floatT* coordinate() {return p.x;}
  floatT x(int i) {return p.x[i];}
  void x(int i, floatT val) {p.x[i]=val;}
  friend bool operator==(pointT a, pointT b) {
    return samePoint<dim>(a.x(), b.x());}
  friend bool operator!=(pointT a, pointT b) {
    return !samePoint<dim>(a.x(), b.x());}
  friend bool operator<(pointT a, pointT b) {
    for(int i=0; i<dim; ++i) {
      if (a[i] != b[i]) {
        if (a[i] < b[i]) return true;
        else return false;
      }}
    return false;
  }
};

template <int dim>
static std::ostream& operator<<(std::ostream& os, const augPoint<dim> v) {
  os << v.p; return os;}

#endif
