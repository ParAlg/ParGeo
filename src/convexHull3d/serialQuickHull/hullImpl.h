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

#include "hullTopology.h"
#include "linkedFacet.h"

#include "parlay/parallel.h"
#include "pargeo/algebra.h"
#include "pargeo/parlayAddon.h"

namespace pargeo::hull3d::serialQuickHull {

  template<class pointT> class hullTopology;

}

template<class pointT>
class pargeo::hull3d::serialQuickHull::hullTopology :
  public pargeo::hull3d::_hullTopology<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>> {

  using facetT = linkedFacet<pointT>;
  using vertexT = pargeo::hull3d::vertex<facetT, pointT>;
  using baseT = _hullTopology<facetT, vertexT>;

 public:

  hullTopology(facetT* f, parlay::sequence<vertexT>& P, vertexT interiorPt):
    baseT() {

    baseT::interiorPt = interiorPt;
    baseT::H = f;

  }

  inline bool keep(facetT* f, vertexT p) {

    if ((f->a - p).dot(f->area) > baseT::eps)
      return f->a != p && f->b != p && f->c != p;
    else
      return false;

  }

  void redistribute(parlay::slice<facetT**, facetT**> facetsBeneath,
		    parlay::slice<facetT**, facetT**> newFacets) {

    baseT::hSize += newFacets.size() - facetsBeneath.size();

    // Redistribute the visible points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numPts(); ++j) { // Point loop
	facetsBeneath[i]->pts(j).seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (keep(newFacets[k], facetsBeneath[i]->pts(j))) {
	    facetsBeneath[i]->pts(j).seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->pts(j));
	    break;
	  }}}}
  }

  facetT* facetWalk() {

    facetT* f = baseT::H;
    size_t fSize = f->numPts();

    auto fVisit = [&](facetT* _f) {return true;};
    auto fDo = [&](facetT* _f) {
      if (_f->numPts() > fSize) {
	fSize = _f->numPts();
	f = _f;
      }
    };
    auto fStop = [&]() { return false; };
    baseT::dfsFacet(f, fVisit, fDo, fStop);
    return f;

  }

  vertexT furthestApex(facetT *f=nullptr) {

    vertexT apex = vertexT();

    auto fVisit = [&](facetT* f) {return true;};
    auto fDo = [&](facetT* f) {
		 if (f->numPts() > 0) apex = f->furthest();
	       };
    auto fStop = [&]() { return !apex.isEmpty(); };
    baseT::dfsFacet(f ? f : baseT::H, fVisit, fDo, fStop);
    return apex;

  }

};

