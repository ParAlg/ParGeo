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

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "convexHull3d/serialQuickHull/hullImpl.h"
#include "convexHull3d/serialQuickHull/linkedFacet.h"
#include "convexHull3d/vertex.h"

// #define HULL_SERIAL_VERBOSE

namespace pargeo {
  namespace hull3d {
    namespace serialQuickHull {

      template<class linkedFacet3d, class vertex3d, class pointT>
      void quickHullSerial(pargeo::hull3d::serialQuickHull::hullTopology<pointT> *context);

    }
  }
}

template<class linkedFacet3d, class vertex3d, class pointT>
void pargeo::hull3d::serialQuickHull::quickHullSerial(pargeo::hull3d::serialQuickHull::hullTopology<pointT> *context) {
  using namespace pargeo;
  using namespace parlay;

#ifdef HULL_SERIAL_VERBOSE
  timer t;
  size_t errors = 0;
  size_t round = 0;
  double apexTime = 0;
  double frontierTime = 0;
  double createTime = 0;
  double splitTime = 0;
#endif

  while (true) {

#ifdef HULL_SERIAL_VERBOSE
    t.start();
#endif

  loopStart:

    linkedFacet3d* f0 = context->hullSize() < 512 ? context->facetWalk() : nullptr;
    vertex3d apex = context->furthestApex(f0);

#ifdef HULL_SERIAL_VERBOSE
    apexTime += t.get_next();
#endif

    if (apex.isEmpty()) break;

#ifdef HULL_SERIAL_VERBOSE
    round ++;
#endif

    auto frontier = context->computeFrontier(apex);
    auto frontierEdges = std::move(std::get<0>(frontier));
    auto facetsBeneath = std::move(std::get<1>(frontier));

#ifdef HULL_SERIAL_VERBOSE
    frontierTime += t.get_next();
#endif

    for(size_t i=0; i<frontierEdges.size(); ++i) {
      auto nv = frontierEdges.at((i+1)%frontierEdges.size());
      auto cv = frontierEdges.at(i);
      if (cv.b != nv.a) {
	apex.seeFacet->clear();

#ifdef HULL_SERIAL_VERBOSE
	errors ++;
#endif

	goto loopStart;
      }
    }

    // Create new facets
    auto newFacets = sequence<linkedFacet3d*>(frontierEdges.size());

    for (size_t i=0; i<frontierEdges.size(); ++i) {
      typename pargeo::hull3d::_hullTopology<linkedFacet3d, vertex3d>::_edge e = frontierEdges.at(i);
      newFacets[i] = new linkedFacet3d(e.a, e.b, apex);
    }

    // Connect new facets
    for (size_t i=0; i<frontierEdges.size(); ++i) {
      context->linkFacet(newFacets[i],
		newFacets[(i+1)%frontierEdges.size()],
		frontierEdges.at(i).ff,
		newFacets[(i-1+frontierEdges.size())%frontierEdges.size()]
		);
    }

    context->setHull(newFacets[0]);

#ifdef HULL_SERIAL_VERBOSE
    createTime += t.get_next();
#endif

    context->redistribute(make_slice(facetsBeneath), make_slice(newFacets));

#ifdef HULL_SERIAL_VERBOSE
    splitTime += t.stop();
#endif

    // Delete existing facets
    for(int j=0; j<facetsBeneath.size(); ++j)
      delete facetsBeneath.at(j);
  }

#ifdef HULL_SERIAL_VERBOSE
  std::cout << "apex-time = " << apexTime << std::endl;
  std::cout << "frontier-time = " << frontierTime << std::endl;
  std::cout << "create-time = " << createTime << std::endl;
  std::cout << "split-time = " << splitTime << std::endl;
  std::cout << "#-rounds = " << round << std::endl;
  std::cout << "#-errors = " << errors << std::endl;
#endif
}
