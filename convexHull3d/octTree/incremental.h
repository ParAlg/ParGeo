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

#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "convexHull3d/hull.h"
#include "gridVertex.h"
#include "gridHull.h"

template<class linkedFacet3d, class vertex3d, class origin3d>
void incrementHull3dSerial(gridHull<linkedFacet3d, vertex3d, origin3d> *context) {
  using namespace pargeo;
  using namespace parlay;

  while (true) {

  loopStart:

    vertex3d apex = context->furthestApex();

    if (apex.isEmpty()) break;

    auto frontier = context->computeFrontier(apex);
    auto frontierEdges = get<0>(frontier);
    auto facetsBeneath = get<1>(frontier);

    for(size_t i=0; i<frontierEdges->size(); ++i) {
      auto nv = frontierEdges->at((i+1)%frontierEdges->size());
      auto cv = frontierEdges->at(i);
      if (cv.b != nv.a) {
	apex.attribute.seeFacet->clear();
	goto loopStart;
      }
    }

    // Create new facets
    auto newFacets = sequence<linkedFacet3d*>(frontierEdges->size());

    for (size_t i=0; i<frontierEdges->size(); ++i) {
      //typename _hullTopology<linkedFacet3d, vertex3d, origin3d>::_edge e = frontierEdges->at(i);
      typename gridHull<linkedFacet3d, vertex3d, origin3d>::_edge e = frontierEdges->at(i);
      newFacets[i] = new linkedFacet3d(e.a, e.b, apex);
    }

    // Connect new facets
    for (size_t i=0; i<frontierEdges->size(); ++i) {
      context->linkFacet(newFacets[i],
		newFacets[(i+1)%frontierEdges->size()],
		frontierEdges->at(i).ff,
		newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		);
    }

    context->setHull(newFacets[0]);

    context->redistribute(facetsBeneath->cut(0, facetsBeneath->size()), make_slice(newFacets));

    // Delete existing facets
    for(int j=0; j<facetsBeneath->size(); ++j)
      delete facetsBeneath->at(j);

    delete frontierEdges;
    delete facetsBeneath;
  }
}
