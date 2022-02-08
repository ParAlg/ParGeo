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

#include "convexHull3d/vertex.h"
#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/serialQuickHull/linkedFacet.h"
#include "convexHull3d/parallelQuickHull/hull.h"
#include "convexHull3d/divideConquer/hull.h"
#include "convexHull3d/sampling/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

// #define HULL_CONCURRENT_VERBOSE

template<class pointT>
using facet1 = pargeo::hull3d::serialQuickHull::linkedFacet<pointT>;

template<class pointT>
using vertex1 = pargeo::hull3d::vertex<facet1<pointT>, pointT>;

template<class pointT>
using facet2 = pargeo::hull3d::parallelQuickHull::linkedFacet<pointT>;

template<class pointT>
using vertex2 = pargeo::hull3d::vertex<facet2<pointT>, pointT>;

template <class pointT>
parlay::sequence<vertex2<pointT>>
concurrentHull(parlay::slice<pointT*, pointT*> P, size_t numProc = 0) {

#ifdef HULL_CONCURRENT_VERBOSE
  pargeo::timer t; t.start();
#endif

  if (!numProc) numProc = parlay::num_workers();

  numProc *= 8;

  size_t blkSize = floor(P.size() / numProc);

  while (blkSize < 10) {
    numProc -= 1;
    blkSize = floor(P.size() / numProc);
  }

#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "#-subproblems = " << numProc << "\n";
#endif

  parlay::sequence<parlay::sequence<vertex2<pointT>>> subHulls(numProc);

  parlay::parallel_for(0, numProc, [&](size_t i) {
				     size_t s = i * blkSize;
				     size_t e = std::min(P.size(), (i+1) * blkSize);

				     parlay::sequence<vertex1<pointT>> S =
				       std::move(pargeo::hull3d::serialQuickHull::computeVertex(P.cut(s, e)));

				     auto cast = [&](size_t i){
						   return vertex2<pointT>(S[i].coords());
						 };

				     subHulls[i] = std::move(parlay::tabulate(S.size(), cast));
				   }, 1);

  parlay::sequence<vertex2<pointT>> uniquePts = parlay::flatten(subHulls);

#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "#-rem-pts = " << uniquePts.size() << "/" << Q.size() << "\n";
  std::cout << "hull-time = " << t.stop() << "\n";
#endif

  // Divide and conquer once seems to be the most efficient for now
  return uniquePts;
}

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::divideConquer::compute(parlay::slice<pointT*, pointT*> P, size_t numProc) {

  if (P.size() < 1000) {
    parlay::sequence<pargeo::hull3d::facet<pointT>> H =
      pargeo::hull3d::serialQuickHull::compute(P);
    return H;
  }

#ifdef HULL_CONCURRENT_VERBOSE
  pargeo::timer t; t.start();
#endif

  // the fraction of points on the hull by sampling
  double frac = pargeo::hull3d::sampling::test<pointT>(P, 0.01);

#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "> sampling-time = " << t.get_next() << "\n";
#endif

  if (frac < 0.8) {

    parlay::sequence<vertex2<pointT>> Q = concurrentHull<pointT>(P, numProc);

#ifdef HULL_CONCURRENT_VERBOSE
    std::cout << "> concurrent-hull-time = " << t.get_next() << "\n";
#endif

    auto H = pargeo::hull3d::parallelQuickHull::compute<pointT>(parlay::make_slice(Q));

#ifdef HULL_CONCURRENT_VERBOSE
    std::cout << "> merge-hull-time = " << t.stop() << "\n";
#endif

    return H;

  } else {

    return pargeo::hull3d::parallelQuickHull::compute<pointT>(parlay::make_slice(P));

  }
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::divideConquer::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, size_t);

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::divideConquer::compute(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, size_t);
