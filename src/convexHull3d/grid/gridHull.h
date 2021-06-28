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

#include "parlay/parallel.h"
#include "pargeo/algebra.h"
#include "pargeo/parlayAddon.h"
#include "convexHull3d/hullTopology.h"

template <class facetT, class vertexT, class originT>
  class gridHull : public _hullTopology<facetT, vertexT, originT> {

  using baseT = _hullTopology<facetT, vertexT, originT>;

 public:

 gridHull(slice<vertexT*, vertexT*> P,
	    originT _origin):
  baseT() {
    baseT::origin = _origin;
    baseT::H = initialize(P);
  }

  facetT* initialize(slice<vertexT*, vertexT*> P) {

    // Maximize triangle area based on fixed xMin and xMax
    size_t X[6];
    auto xx = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[0]<j[0];});
    X[0] = xx.first - &P[0]; X[1] = xx.second - &P[0];
    auto yy = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[1]<j[1];});
    X[2] = yy.first - &P[0]; X[3] = yy.second - &P[0];
    auto zz = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[2]<j[2];});
    X[4] = zz.first - &P[0]; X[5] = zz.second - &P[0];

    size_t xMin, xMax;
    if (P[X[1]][0]-P[X[0]][0] > P[X[3]][1]-P[X[2]][1] && P[X[1]][0]-P[X[0]][0] > P[X[5]][2]-P[X[4]][2]) {
      xMin = X[0]; xMax = X[1];
    } else if (P[X[3]][1]-P[X[2]][1] > P[X[1]][0]-P[X[0]][0] && P[X[3]][1]-P[X[2]][1] > P[X[5]][2]-P[X[4]][2]) {
      xMin = X[2]; xMax = X[3];
    } else {
      xMin = X[4]; xMax = X[5];
    }

    vertexT x1 = P[xMin];
    vertexT x2 = P[xMax];

    auto y = max_element_serial(P, [&](vertexT i, vertexT j) {
	return crossProduct3d(x1-i, x2-i).length() <
	  crossProduct3d(x1-j, x2-j).length();
      });
    size_t yApex = y - &P[0];
    vertexT y1 = P[yApex];

    // Maximize simplex volume
    vertexT area = crossProduct3d(x1-y1, x2-y1);
    auto z = max_element(P, [&](vertexT i, vertexT j) {
	return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
      });
    size_t zApex = z - &P[0];

    size_t c1 = xMin;
    size_t c2 = xMax;
    size_t c3 = yApex;
    size_t c4 = zApex;

    baseT::hSize = 4;

    baseT::origin.setOrigin((P[c1] + P[c2] + P[c3] + P[c4])/4);

    // Initialize points with visible facet link
    auto Q = typename facetT::seqT(P.size());

    for (size_t i=0; i<P.size(); ++i)
      Q[i] = P[i] - baseT::origin.get(); // translation

    // Make initial facets
    auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
    auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
    auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
    auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

    baseT::linkFacet(f0, f1, f2, f3);
    baseT::linkFacet(f1, f0, f2, f3);
    baseT::linkFacet(f2, f1, f0, f3);
    baseT::linkFacet(f3, f1, f2, f0);

    for(size_t i=0; i<Q.size(); i++) {
      if (baseT::origin.visible(f0, Q[i])) {
	Q[i].attribute.seeFacet = f0;
	f0->push_visible(Q[i]);
      } else if (baseT::origin.visible(f1, Q[i])) {
	Q[i].attribute.seeFacet = f1;
	f1->push_visible(Q[i]);
      } else if (baseT::origin.visible(f2, Q[i])) {
	Q[i].attribute.seeFacet = f2;
	f2->push_visible(Q[i]);
      } else if (baseT::origin.visible(f3, Q[i])) {
	Q[i].attribute.seeFacet = f3;
	f3->push_visible(Q[i]);
      } else {
	Q[i].attribute.seeFacet = nullptr;
      }
    }

    for(size_t i=0; i<Q.size(); i++) {
      if (baseT::origin.keep(f0, Q[i])) f0->push_keep(Q[i]);
      if (baseT::origin.keep(f1, Q[i])) f1->push_keep(Q[i]);
      if (baseT::origin.keep(f2, Q[i])) f2->push_keep(Q[i]);
      if (baseT::origin.keep(f3, Q[i])) f3->push_keep(Q[i]);
    }

    return f0;
  }

  void redistribute(slice<facetT**, facetT**> facetsBeneath,
		    slice<facetT**, facetT**> newFacets) {

    baseT::hSize += newFacets.size() - facetsBeneath.size();

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numVisiblePts(); ++j) { // Point loop
	facetsBeneath[i]->visiblePts(j).attribute.seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (baseT::origin.visible(newFacets[k], facetsBeneath[i]->visiblePts(j))) {
	    facetsBeneath[i]->visiblePts(j).attribute.seeFacet = newFacets[k];
	    newFacets[k]->push_visible(facetsBeneath[i]->visiblePts(j));
	    break;
	  }}}}

    sequence<vertexT> intersection;
    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numKeepPts(); ++j) { // Point loop
	intersection.push_back(facetsBeneath[i]->keepPts(j));
      }}
    parlay::sort_inplace(intersection); // todo make serial
    intersection = unique(intersection);
    for (auto p: intersection) {
      p.attribute.seeFacet = nullptr;
      for (int k=0; k<nnf; ++k) { // New facet loop
	if (baseT::origin.keep(newFacets[k], p))
	  newFacets[k]->push_keep(p);
      }}

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
