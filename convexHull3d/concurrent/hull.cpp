#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "pargeo/zorderSort.h"
#include "convexHull3d/pointVertex.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/incremental.h"
#include "convexHull3d/hull.h"

using namespace pargeo;

template <typename ptOut>
parlay::sequence<ptOut> concurrentHull(parlay::sequence<pointVertex> &Q, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using pt = pointVertex;
#ifndef SILENT
  cout << "------------" << endl;
  cout << "input-size = " << Q.size() << endl;
  timer t; t.start();
#endif

  //sort_inplace(Q, [&](pt i, pt j){return i[0] < j[0];});
  //zorderSortInPlace3d(Q);

  if (!numProc) numProc = num_workers();

  numProc *= 8;

  size_t blkSize = floor(Q.size() / numProc);

  while (blkSize < 100) {
    numProc -= 1;
    blkSize = floor(Q.size() / numProc);
  }
#ifndef SILENT
  cout << "#-subproblems = " << numProc << endl;
#endif
  sequence<sequence<ptOut>> subHulls(numProc);

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = min(Q.size(), (i+1) * blkSize);
			     auto origin = pointOrigin();
			     auto linkedHull = new _hull<linkedFacet3d<pt>, pt, pointOrigin>(Q.cut(s, e), origin, true);
			     incrementHull3dSerial<linkedFacet3d<pt>, pt, pointOrigin>(linkedHull);
			     subHulls[i] = linkedHull->getHullVertices<ptOut>();
			     //subHulls[i] = hull3dInternalSerial(Q.cut(s, e)); // todo compilation issue
			   }, 1);

  sequence<ptOut> uniquePts = parlay::flatten(subHulls);
#ifndef SILENT
  cout << "output-size = " << uniquePts.size() << endl;
  cout << "hull-time = " << t.stop() << endl;
#endif

  // Divide and conquer once seems to be the most efficient for now
  return uniquePts;

  // note can consider recursively calling concurrent hull, did not find fast before
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3dConcurrent(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using pt = pointVertex;

  size_t n = P.size();

#ifndef SILENT
  cout << "#-points = " << n << endl;
  timer t; t.start();
#endif

  sequence<pt> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = pt(P[i].coords());});

#ifndef SILENT
  cout << "> concurrent-hull-time = " << t.get_next() << endl;
#endif

  sequence<pt> Q2 = concurrentHull<pt>(Q, numProc);
  auto origin = pointOrigin();
  auto finalLinkedHull = new _hull<linkedFacet3d<pt>, pt, pointOrigin>(make_slice(Q2), origin, false);
  incrementHull3dSerial<linkedFacet3d<pt>, pt, pointOrigin>(finalLinkedHull); // todo parallelize
  //incrementHull3d<linkedFacet3d<pt>, pt, pointOrigin>(finalLinkedHull, numProc);
  auto out = sequence<facet3d<pargeo::fpoint<3>>>();
  finalLinkedHull->getHull<pargeo::fpoint<3>>(out);

#ifndef SILENT
  cout << "> merge-hull-time = " << t.stop() << endl;
#endif

#ifndef SILENT
  cout << "hull-size = " << out.size() << endl;
#endif

  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3dConcurrent(parlay::sequence<pargeo::fpoint<3>> &, size_t);
