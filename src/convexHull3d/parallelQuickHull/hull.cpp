
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

#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "convexHull3d/initialize.h"
#include "convexHull3d/bruteforce/hull.h"

#include "convexHull3d/parallelQuickHull/hull.h"
#include "convexHull3d/parallelQuickHull/hullImpl.h"
#include "convexHull3d/parallelQuickHull/quickHull.h"
#include "convexHull3d/parallelQuickHull/linkedFacet.h"

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<pointT*, pointT*> P, size_t numProc, bool randomized) {

  hullTopology<pointT>* linkedHull = pargeo::hull3d::initParallel<hullTopology<pointT>, linkedFacet<pointT>, pointT>(P);

  quickHullParallel<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>>(linkedHull, numProc, randomized);

  auto out = parlay::sequence<pargeo::hull3d::facet<pointT>>();
  linkedHull->getFacet(out);

  delete linkedHull;
  return out;
}

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<
					   pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pointT>, pointT>*,
					   pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pointT>, pointT>*
					   > P, size_t numProc) {

  hullTopology<pointT>* linkedHull = pargeo::hull3d::initParallel<hullTopology<pointT>, linkedFacet<pointT>, pointT>(P);

  quickHullParallel<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>>(linkedHull, numProc);

  auto out = parlay::sequence<pargeo::hull3d::facet<pointT>>();
  linkedHull->getFacet(out);

  delete linkedHull;
  return out;
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<
	pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pargeo::fpoint<3>>, pargeo::fpoint<3>>*,
	pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pargeo::fpoint<3>>, pargeo::fpoint<3>>*
	>, size_t = 0);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<
	pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pargeo::point<3>>, pargeo::point<3>>*,
	pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pargeo::point<3>>, pargeo::point<3>>*
	>, size_t = 0);

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::parallelQuickHull::compute(hullTopology<pointT>* linkedHull, size_t numProc) {

  quickHullParallel<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>>(linkedHull, numProc);

  auto out = parlay::sequence<pargeo::hull3d::facet<pointT>>();
  linkedHull->getFacet(out);

  return out;
}

// parlay::sequence<pargeo::facet<pargeo::fpoint<3>>>
// pargeo::hullInternal::hull3dParallelInternal(parlay::slice<
// 			       pargeo::hullInternal::vertex*,
// 			       pargeo::hullInternal::vertex*
// 			       > Q,
// 			       size_t numProc) {
//   using namespace std;
//   using namespace parlay;
//   using floatT = pargeo::fpoint<3>::floatT;
//   using pointT = pargeo::fpoint<3>;
//   using facetT = facet<pargeo::fpoint<3>>;
//   using vertexT = pargeo::hullInternal::vertex;

//   // Create an initial simplex
//   auto origin = pointOrigin();

//   auto linkedHull = new hullTopology<linkedFacet<vertexT>, vertexT, pointOrigin>(Q, origin);

//   incrementHull3d<linkedFacet<vertexT>, vertexT, pointOrigin>(linkedHull, numProc);

//   // getHull will undo the translation of linkedHull
//   auto out = sequence<facetT>();
//   linkedHull->getHull<pointT>(out);

//   delete linkedHull;
//   return out;
// }

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, size_t, bool);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::parallelQuickHull::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, size_t, bool);
