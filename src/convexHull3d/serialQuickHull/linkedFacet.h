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

#include <vector>
#include "convexHull3d/vertex.h"
#include "pargeo/algebra.h"

namespace pargeo::hull3d::serialQuickHull {

  template<class pointT> class linkedFacet;

}

template<class pointT>
class pargeo::hull3d::serialQuickHull::linkedFacet {

public:

  using vertexT = pargeo::hull3d::vertex<linkedFacet, pointT>;
  using floatT = typename vertexT::floatT;

  vertexT a, b, c;

  linkedFacet *abFacet;

  linkedFacet *bcFacet;

  linkedFacet *caFacet;

  vertexT area;

  std::vector<vertexT> seeList;

  size_t numPts() { return seeList.size(); }

  vertexT& pts(size_t i) { return seeList.at(i); }

  void clear() { seeList.clear(); }

  void push_back(vertexT v) { seeList.push_back(v); }

  vertexT furthest() {
    auto apex = vertexT();
    typename vertexT::floatT m = pointT::eps;
    for (size_t i=0; i<numPts(); ++i) {
      auto m2 = (a-pts(i)).dot(area);
      if (m2 > m) {
	m = m2;
	apex = pts(i);
      }
    }
    return apex;
  }

  vertexT furthestSample(size_t s) {
    auto apex = vertexT();
    typename vertexT::floatT m = pointT::eps;
    for (size_t i=0; i<std::min(numPts(), s); ++i) {
      size_t ii = parlay::hash64(i) % numPts();
      auto m2 = (a-pts(ii)).dot(area);
      if (m2 > m) {
	m = m2;
	apex = pts(ii);
      }
    }
    return apex;
  }

  linkedFacet(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > pointT::eps)
      std::swap(b, c);

    seeList = std::vector<vertexT>();

    area = crossProduct3d(b-a, c-a);
  }

};
