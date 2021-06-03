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

#include "incremental.h"
#include "gridVertex.h"
#include "gridHull.h"
#include "octTree.h"

using namespace pargeo;
using namespace parlay;

parlay::sequence<pargeo::fpoint<3>>
hull3dGridInternal(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*> P,
		   size_t s = 4) {
  using namespace std;
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pargeo::fpoint<3>>;

  auto levels = sequence<size_t>();
  levels.push_back(s);
  //levels.push_back(15);
  //levels.push_back(4);
  //levels.push_back(2);
  auto tree = octTree<gridVertex>(P, levels);

  size_t L = tree.numLevels();

  sequence<gridVertex> Q = tree.levelSample(0); // Get the coarsest level sample

  sequence<size_t> mask;

  for (size_t l = 0; l < tree.numLevels(); ++l) {

    if (l > 0) Q = tree.nextLevelSample(l-1, mask);

    auto origin = gridOrigin();

    origin.setMin(tree.coordinateMin());

    origin.setBoxSize(tree.boxSize(l));

    auto linkedHull = new gridHull<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(make_slice(Q), origin);

    origin = linkedHull->getOrigin();

    incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(linkedHull);

    auto pts = linkedHull->getHullPts();

    mask = sequence<size_t>(tree.levelSize(l) + 1, 0);
    parallel_for(0, pts.size(), [&](size_t i) {
		   mask[pts[i].attribute.i] = 1;
		 });

    delete linkedHull;
  }

  return tree.getPoints<pointT>(tree.numLevels()-1, mask);
}

parlay::sequence<facet3d<pargeo::fpoint<3>>>
pargeo::hull3dGrid(parlay::sequence<pargeo::fpoint<3>> &P, size_t s = 4, bool write = false) {
  using namespace std;
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pargeo::fpoint<3>>;

  pargeo::timer t;
  t.start();

  auto finalQ = hull3dGridInternal(make_slice(P), s);

  cout << " " << finalQ.size() << " / " << P.size() << endl;
  cout << "grid-par-time = " << t.get_next() << endl;

  // todo change merge hull method here
  auto result = hull3dConcurrent(finalQ);
  cout << "merge-hull-time = " << t.stop() << endl;

  cout << "hull-size = " << result.size() << endl;
  return result;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>>
pargeo::hull3dGridConcurrent(parlay::sequence<pargeo::fpoint<3>> &P, size_t s, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pargeo::fpoint<3>>;

  pargeo::timer t;

  t.start();

  if (!numProc) numProc = num_workers();

  numProc *= 8;

  size_t blkSize = floor(P.size() / numProc);

  while (blkSize < 100) {
    numProc -= 1;
    blkSize = floor(P.size() / numProc);
  }

  cout << "#-subproblems = " << numProc << endl;

  sequence<sequence<pointT>> subHulls(numProc);

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = min(P.size(), (i+1) * blkSize);
			     cout << s << " -- " << e << endl;
			     subHulls[i] = hull3dGridInternal(P.cut(s, e));
			   }, 1);

  sequence<pointT> finalQ = parlay::flatten(subHulls);

  cout << " " << finalQ.size() << " / " << P.size() << endl;
  cout << "grid-par-time = " << t.get_next() << endl;

  // todo change merge method here
  auto result = hull3dConcurrent(finalQ);

  cout << "merge-hull-time = " << t.stop() << endl;

  cout << "hull-size = " << result.size() << endl;
  return result;
}
