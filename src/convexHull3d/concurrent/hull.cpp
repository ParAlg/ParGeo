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

#include "convexHull3d/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace pargeo;

// #define HULL_CONCURRENT_VERBOSE

template <typename ptOut>
parlay::sequence<ptOut>
concurrentHull(parlay::sequence<vertex> &Q, size_t numProc) {
  using namespace parlay;
  using pt = vertex;

#ifdef HULL_CONCURRENT_VERBOSE
  timer t; t.start();
#endif

  if (!numProc) numProc = num_workers();

  numProc *= 8;

  size_t blkSize = floor(Q.size() / numProc);

  while (blkSize < 10) {
    numProc -= 1;
    blkSize = floor(Q.size() / numProc);
  }

#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "#-subproblems = " << numProc << "\n";
#endif

  sequence<sequence<ptOut>> subHulls(numProc);

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = std::min(Q.size(), (i+1) * blkSize);
			     subHulls[i] =
			       std::move(hullInternal::hull3dSerialInternal(Q.cut(s, e)));
			   }, 1);

  sequence<ptOut> uniquePts = parlay::flatten(subHulls);
#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "#-rem-pts = " << uniquePts.size() << "/" << Q.size() << "\n";
  std::cout << "hull-time = " << t.stop() << "\n";
#endif

  // Divide and conquer once seems to be the most efficient for now
  return uniquePts;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>>
pargeo::hull3dConcurrent(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using pt = vertex;

  if (P.size() < 1000) return hull3dSerial(P);

  size_t n = P.size();

#ifdef HULL_CONCURRENT_VERBOSE
  timer t; t.start();
#endif

  // the fraction of points on the hull by sampling
  double frac = testHull(P, 0.01);

#ifdef HULL_CONCURRENT_VERBOSE
  std::cout << "> sampling-time = " << t.get_next() << "\n";
#endif

  if (frac < 0.8) {

    sequence<pt> Q(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				Q[i] = pt(P[i].coords());});

    sequence<pt> Q2 = concurrentHull<pt>(Q, numProc);

#ifdef HULL_CONCURRENT_VERBOSE
    std::cout << "> concurrent-hull-time = " << t.get_next() << "\n";
#endif

    return hull3dSerialInternal(make_slice(Q2));
    // return hull3dIncrementalInternal(make_slice(Q2));

#ifdef HULL_CONCURRENT_VERBOSE
    std::cout << "> merge-hull-time = " << t.stop() << "\n";
#endif

  } else {

    return hull3dSerial(P);

  }
}
