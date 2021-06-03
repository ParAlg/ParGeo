// This code is part of the Pargeo Library
// Copyright (c) 2021 Yiqiu Wang and the Pargeo Team
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
#include <iostream>
#include <fstream>
#include <bitset>
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

class gridAtt3d;

using gridVertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, gridAtt3d>;

template <class vertexT> struct linkedFacet3d;

struct gridAtt3d {
  gridAtt3d() {}

  // ------------ oct tree data fields

  size_t id; // multi-level grid id

  size_t i; // index in the level

  // ------------ convex hull data fields

  linkedFacet3d<gridVertex> *seeFacet;
};

static std::ostream& operator<<(std::ostream& os, const gridVertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

class gridOrigin;

template <class vertexT>
struct linkedFacet3d {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

  //#ifdef SERIAL_HULL
  //typedef std::vector<vertexT> seqT;
  //#else
  typedef parlay::sequence<vertexT> seqT;
  //#endif

  vertexT a, b, c;
  linkedFacet3d *abFacet;
  linkedFacet3d *bcFacet;
  linkedFacet3d *caFacet;
  vertexT area;

  bool isAdjacent(linkedFacet3d* f2) {
    vertexT V[3]; V[0] = a; V[1] = b; V[2] = c;
    for (int i=0; i<3; ++i) {
      auto v1 = V[i];
      auto v2 = V[(i+1)%3];
      if ( (f2->a == v1 && f2->b == v2) || (f2->a == v2 && f2->b == v1) ) return true;
      if ( (f2->b == v1 && f2->c == v2) || (f2->b == v2 && f2->c == v1) ) return true;
      if ( (f2->c == v1 && f2->a == v2) || (f2->c == v2 && f2->a == v1) ) return true;
    }
    return false;
  }

  seqT *seeList;

  seqT *keepList;

  void reassignVisible(seqT *_seeList) {
    delete seeList;
    seeList = _seeList;
  }

  void reassignKeep(seqT *_keepList) {
    delete keepList;
    keepList = _keepList;
  }

  void push_back(vertexT v, gridOrigin* o) {
    if (o->visible(this, v) && a != v && b != v && c != v)
      seeList->push_back(v);
    else
      keepList->push_back(v);
  }

  void push_keep(vertexT v) {
    keepList->push_back(v);
  }

  void push_visible(vertexT v) {
    if (a != v && b != v && c != v)
      seeList->push_back(v);
  }

  // Accesses the visible samples

  size_t numVisiblePts() { return seeList->size(); }

  inline vertexT& visiblePts(size_t i) { return seeList->at(i); }

  // Accesses the boxes

  size_t numKeepPts() { return keepList->size(); }

  inline vertexT& keepPts(size_t i) { return keepList->at(i); }

  // Access both the visible samples and intersecting boxes

  inline size_t numPts() { return seeList->size() + keepList->size(); }

  vertexT& pts(size_t i) {
    if (i < seeList->size()) {
      return seeList->at(i);
    } else {
      return keepList->at(i - seeList->size());
    }
  }

  void clear() {
    seeList->clear();
    keepList->clear();
  }

  vertexT furthest() {
    auto apex = vertexT();
    typename vertexT::floatT m = numericKnob;

    for (size_t i=0; i<numVisiblePts(); ++i) {
      auto m2 = pargeo::signedVolumeX6(a, b, c, visiblePts(i));
      if (m2 > m) {
	m = m2;
	apex = visiblePts(i);
      }
    }
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > numericKnob)
      std::swap(b, c);

    seeList = new seqT();
    keepList = new seqT();

    area = crossProduct3d(b-a, c-a);
  }

  ~linkedFacet3d() {
    delete seeList;
    delete keepList;
  }
};

template<class vertexT>
static std::ostream& operator<<(std::ostream& os, const linkedFacet3d<vertexT>& v) {
#ifdef VERBOSE
  os << "(" << v.a.attribute.i << "," << v.b.attribute.i << "," << v.c.attribute.i << ")";
#else
  os << "(" << v.a << "," << v.b << "," << v.c << ")";
#endif
  return os;
}

class gridOrigin {
  using facetT = linkedFacet3d<gridVertex>;
  using vertexT = gridVertex;
  using pt = pargeo::fpoint<3>;
  using floatT = typename pt::floatT;

  static constexpr floatT numericKnob = 1e-5;

  vertexT o; // translation origin

  vertexT pMin; // translated pmin of the grid

  floatT bSize; // gridSize of the level being processed

public:
  void setBoxSize(floatT _boxSize) {
    bSize = _boxSize;
  }

  template<class pt>
  void setMin(pt _pMin) {
    pMin = vertexT(_pMin.coords());
  }

  template<class pt>
  void setOrigin(pt _o) {
    o = vertexT(_o.coords());
    pMin = pMin - o; // translate the pmin by origin (same as other points in the hull)
  }

  inline vertexT get() {return o;};

  gridOrigin() {}

  template<class pt>
  gridOrigin(pt _o) {
    setOrigin(_o);
  }

  // Point visibility test
  inline bool visible(facetT* f, vertexT p) {
    return pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob;
  }

  vertexT getCorner(vertexT p, vertexT myMin) {
    // todo compute from point grid id
    vertexT vt0 = p - myMin;
    vt0[0] = floor( vt0[0] / bSize ) * bSize;
    vt0[1] = floor( vt0[1] / bSize ) * bSize;
    vt0[2] = floor( vt0[2] / bSize ) * bSize;
    return vt0 + myMin;
  }

  // // Input p needs to be untranslated
  // inline void writeCorners(vertexT p, std::ofstream& os) {
  //   if (bSize < 0) {
  //     os << p[0] << " " << p[1] << " " << p[2] << " " << 0 << std::endl;
  //   } else {
  //     vertexT vt0 = getCorner(p, pMin + o);
  //     for (size_t i=0; i<8; ++i) {
  // 	auto vt = vt0;
  // 	vt[0] += (i & 1) * bSize;
  // 	vt[1] += ((i>>1) & 1) * bSize;
  // 	vt[2] += ((i>>2) & 1) * bSize;
  // 	os << vt[0] << " " << vt[1] << " " << vt[2] << " " << i << std::endl;
  //     }
  //   }
  // }

  inline bool anyCornerVisible(facetT* f, vertexT p) {
    if (bSize < 0) {
      return visible(f, p);
    } else {
      vertexT vt0 = getCorner(p, pMin);
      for (size_t i=0; i<8; ++i) {
	vertexT vt = vt0;
	vt[0] += (i & 1) * bSize;
	vt[1] += ((i>>1) & 1) * bSize;
	vt[2] += ((i>>2) & 1) * bSize;

	if (visible(f, vt))
	  return true;
      }
      return false;
    }
  }

  // This function returns true of any of the box corner of p is visible
  inline bool keep(facetT* f, vertexT p) {
    return anyCornerVisible(f, p);
  }
};
