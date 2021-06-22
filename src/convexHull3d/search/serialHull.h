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

#include "convexHull3d/hullTopology.h"

#include <limits>
#include "parlay/parallel.h"
#include "pargeo/algebra.h"
#include "pargeo/parlayAddon.h"
#include "pargeo/kdTree.h"

template <class facetT, class vertexT, class originT>
class serialHull : public _hullTopology<facetT, vertexT, originT> {

  using baseT = _hullTopology<facetT, vertexT, originT>;

  using floatT = float;

  using nodeT = pargeo::kdNode<3, vertexT>;

 public:
  nodeT* tree;

  serialHull(slice<vertexT*, vertexT*> P,
	     originT _origin):
    baseT() {
    baseT::origin = _origin;
    baseT::H = initialize(P);
  }

  sequence<vertexT> Q;

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
    Q = typename facetT::seqT(P.size());
    //auto Q = typename facetT::seqT(P.size());

    for (size_t i=0; i<P.size(); ++i)
      Q[i] = P[i] - baseT::origin.get(); // translation

    tree = pargeo::buildKdt<3, vertexT>(Q, true, false);

    // Make initial facets
    auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
    auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
    auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
    auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

    baseT::linkFacet(f0, f1, f2, f3);
    baseT::linkFacet(f1, f0, f2, f3);
    baseT::linkFacet(f2, f1, f0, f3);
    baseT::linkFacet(f3, f1, f2, f0);

    // for(size_t i=0; i<Q.size(); i++) {
    //   if (baseT::origin.keep(f0, Q[i])) {
    // 	Q[i].attribute.seeFacet = f0;
    // 	f0->push_back(Q[i]);
    //   } else if (baseT::origin.keep(f1, Q[i])) {
    // 	Q[i].attribute.seeFacet = f1;
    // 	f1->push_back(Q[i]);
    //   } else if (baseT::origin.keep(f2, Q[i])) {
    // 	Q[i].attribute.seeFacet = f2;
    // 	f2->push_back(Q[i]);
    //   } else if (baseT::origin.keep(f3, Q[i])) {
    // 	Q[i].attribute.seeFacet = f3;
    // 	f3->push_back(Q[i]);
    //   } else {
    // 	Q[i].attribute.seeFacet = nullptr;
    //   }
    // }

    return f0;
  }

  void redistribute(slice<facetT**, facetT**> facetsBeneath,
		    slice<facetT**, facetT**> newFacets) {

    baseT::hSize += newFacets.size() - facetsBeneath.size();
    return;

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numPts(); ++j) { // Point loop
	facetsBeneath[i]->pts(j).attribute.seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (baseT::origin.keep(newFacets[k], facetsBeneath[i]->pts(j))) {
	    facetsBeneath[i]->pts(j).attribute.seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->pts(j));
	    break;
	  }}}}
  }

  std::tuple<floatT, floatT> nodeVol(nodeT* r, facetT* f) {
    floatT vMin = std::numeric_limits<floatT>::max();
    floatT vMax = -1;
    floatT vol;
    vertexT c;

    c[0] = r->getMin(0); c[1] = r->getMin(1); c[2] = r->getMin(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMax(0); c[1] = r->getMin(1); c[2] = r->getMin(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMin(0); c[1] = r->getMax(1); c[2] = r->getMin(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMin(0); c[1] = r->getMin(1); c[2] = r->getMax(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMax(0); c[1] = r->getMax(1); c[2] = r->getMin(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMin(0); c[1] = r->getMax(1); c[2] = r->getMax(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMax(0); c[1] = r->getMin(1); c[2] = r->getMax(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    c[0] = r->getMax(0); c[1] = r->getMax(1); c[2] = r->getMax(2);
    vol = f->getVolume(c);
    vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
    return std::tuple(vMin, vMax);
  }

  /*
    - r keeps the furthest distance found for any visible points
   */
  void searchHelper(nodeT* r, facetT* f, float& d, vertexT& apex, size_t& counter) {

    auto vols = nodeVol(r, f);
    floatT volMin = std::get<0>(vols);
    floatT volMax = std::get<1>(vols);

    if ((volMin <= 0 && volMax <= 0) || volMax <= d) {
      // node invisible
      // std::cout << "invisible\n";
      return;
    } else if (volMin <= 0 && volMax > 0) {
      // node intersect
      //std::cout << "intersect\n";
      if (!r->isLeaf()) {
	searchHelper(r->L(), f, d, apex, counter);
	searchHelper(r->R(), f, d, apex, counter);
      } else {
	for (size_t i = 0; i < r->size(); ++ i) {
	  counter ++;
	  auto pVol = f->getVolume(*r->at(i));
	  if (pVol > d) {
	    d = pVol;
	    apex = *r->at(i);
	    //std::cout << " " << *r->at(i) << ": " << pVol << "\n";
	  }
	}
      }
    } else {
      // node visible
      //std::cout << "visible\n";
      for (size_t i = 0; i < r->size(); ++ i) {
	counter ++;
        auto pVol = f->getVolume(*r->at(i));
	if (pVol > d) {
	  d = pVol;
	  apex = *r->at(i);
	}
      }
    }

  }

  vertexT searchVisible(facetT* f) {
    floatT d = 1e-5;
    vertexT apex;
    size_t counter = 0;
    apex.setEmpty();
    searchHelper(tree, f, d, apex, counter);
    if (f->a == apex ||
	f->b == apex ||
	f->c == apex) {
      return vertexT();
    }
    //std::cout << "search counter = " << counter << "\n";
    return apex;
  }

  vertexT searchBrute(facetT* f) {
    size_t counter = 0;
    floatT d = 1e-5;
    vertexT apex;
    apex.setEmpty();
    for (auto p: Q) {
      auto pVol = f->getVolume(p);
      if (pVol > 0) counter ++;
      if (pVol > d &&
	  f->a != p &&
	  f->b != p &&
	  f->c != p) {
	//std::cout << "assign " << pVol << " from " << p << " \n";
	d = pVol;
	apex = p;
      }
    }
    std::cout << "brute counter = " << counter << "\n";
    return apex;
  }

  vertexT furthestApex(facetT *f=nullptr) {
    vertexT apex = vertexT();

    //std::cout << "compute apex\n";
    auto fVisit = [&](facetT* f) {return true;};
    auto fDo = [&](facetT* f) {
		 // todo call a tree query here instead
		 // if (f->numPts() > 0) {
		 //   apex = f->furthest();
		 //   std::cout << apex << ": " << f->getVolume(apex) << "\n";
		 //   auto apex2 = searchBrute(f);
		 //   std::cout << apex2 << ": " << f->getVolume(apex2) << "\n";
		 //   auto apex3 = searchVisible(f);
		 //   std::cout << apex3 << ": " << f->getVolume(apex3) << "\n";
		 //   abort();
		 // }
		 //std::cout << "search " << f << "\n";
		 // apex = searchBrute(f);
		 // auto apex2 = searchVisible(f);
		 // if (apex != apex2) {
		 //   std::cout << "search diff: \n";
		 //   std::cout << *f << "\n";
		 //   std::cout << "brute = " << apex << ": " << f->getVolume(apex) << "\n";
		 //   std::cout << "tree = " << apex2 << ": " << f->getVolume(apex2) << "\n";
		 //   abort();
		 // }
		 if (f->hasVisible) {
		   //auto apexBrute = searchBrute(f);
		   apex = searchVisible(f);
		   // std::cout << "\n";
		   if (!apex.isEmpty()) {
		     //std::cout << apex << ": " << f->getVolume(apex) << "\n";
		     apex.attribute.seeFacet = f;
		   } else {
		     f->hasVisible = false;
		   }
		 }
	       };
    auto fStop = [&]() { return !apex.isEmpty(); };

    baseT::dfsFacet(f ? f : baseT::H, fVisit, fDo, fStop);
    return apex;
  }

  // Also checks the hull integrity
  void printFacets(facetT* start=nullptr) {

    std::cout << "hull-facets = ";
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 cout << f << " ";
	       };
    auto fStop = [&]() { return false;};

    baseT::dfsFacet(baseT::H, fVisit, fDo, fStop);
    std::cout << "\n";
  }

};
