// This code is part of the Pargeo Library
// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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

#pragma once

//#include <vector>
#include "convexHull3d/searchHull.h"
#include "parlay/sequence.h"
#include "pargeo/algebra.h"

class vertexAtt;

using vertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt>;

template <class vertexT> struct linkedFacet3d;

// should probably template the facet
class vertexAtt {
public:
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

// #ifdef VERBOSE
//   size_t i;
// #endif
  linkedFacet3d<vertex> *seeFacet;
  vertexAtt() {}
};

static std::ostream& operator<<(std::ostream& os, const vertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template <class vertexT>
struct linkedFacet3d {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;
  typedef parlay::sequence<vertexT> seqT;

  vertexT a, b, c;

  linkedFacet3d *abFacet;

  linkedFacet3d *bcFacet;

  linkedFacet3d *caFacet;

  vertexT area;

  // seqT *seeList;
  seqT seeList; // todo can remove
  bool hasVisible;

  size_t numPts() {
    return seeList.size(); }

  vertexT& pts(size_t i) {
    return seeList.at(i); }

  void clear() {
    seeList.clear(); }

  void push_back(vertexT v) {
    seeList.push_back(v); }

  float getVolume(vertexT p) { // todo type
    return (a-p).dot(area);
  }

  vertexT furthest() {
    auto apex = vertexT();
    typename vertexT::floatT m = numericKnob;
    for (size_t i=0; i<numPts(); ++i) {
      //auto m2 = (a-pts(i)).dot(crossProduct3d(b-a, c-a));
      auto m2 = (a-pts(i)).dot(area);
      if (m2 > m) {
	m = m2;
	apex = pts(i);
      }
    }
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > numericKnob)
      std::swap(b, c);

    hasVisible = true;
    //seeList = new seqT();
    seeList = seqT(); // todo can remove

    area = crossProduct3d(b-a, c-a);
  }

  ~linkedFacet3d() {
    //delete seeList;
  }

};

template<class vertexT>
static std::ostream& operator<<(std::ostream& os, const linkedFacet3d<vertexT>& v) {
#ifdef HULL_SERIAL_VERBOSE
  os << "(" << v.a.attribute.i << "," << v.b.attribute.i << "," << v.c.attribute.i << ")";
#else
  os << "(" << v.a << "," << v.b << "," << v.c << ")";
#endif
  return os;
}

class pointOrigin {
  using facetT = linkedFacet3d<vertex>;
  using vertexT = vertex;

  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

  vertexT o;

public:

  void setOrigin(vertexT _o) { o = _o; }

  inline vertexT get() {return o;};

  pointOrigin() {}

  template<class pt>
  pointOrigin(pt _o) {
    o = vertexT(_o.coords());
  }

  inline bool visible(facetT* f, vertexT p) {
    return (f->a - p).dot(f->area) > numericKnob &&
      f->a != p && f->b != p && f->c != p;
  }

  inline bool keep(facetT* f, vertexT p) {
    if ((f->a - p).dot(f->area) > numericKnob)
      return f->a != p && f->b != p && f->c != p;
    else
      return false;
  }

};
