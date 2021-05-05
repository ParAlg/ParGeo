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

parlay::sequence<pointVertex> concurrentHull(parlay::sequence<pointVertex> &Q, size_t numProc = 0) {
  using namespace std;
  using namespace parlay;
  using pt = pointVertex;

  cout << "------------" << endl;
  cout << "input-size = " << Q.size() << endl;

  timer t; t.start();

  //sort_inplace(Q, [&](pt i, pt j){return i[0] < j[0];});
  //zorderSortInPlace3d(Q);

  if (!numProc) numProc = num_workers();

  numProc *= 8;

  size_t blkSize = floor(Q.size() / numProc);

  while (blkSize < 100) {
    numProc -= 1;
    blkSize = floor(Q.size() / numProc);
  }

  cout << "#-subproblems = " << numProc << endl;

  sequence<sequence<pt>> subHulls(numProc);

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = min(Q.size(), (i+1) * blkSize);
			     auto origin = pointOrigin();
			     auto linkedHull = new _hull<linkedFacet3d<pt>, pt, pointOrigin>(Q.cut(s, e), origin, true);
			     incrementHull3dSerial<linkedFacet3d<pt>, pt, pointOrigin>(linkedHull);
			     subHulls[i] = linkedHull->getHullVertices<pt>();
			   }, 1);

  sequence<pt> uniquePts = parlay::flatten(subHulls);

  cout << "output-size = " << uniquePts.size() << endl;
  cout << "hull-time = " << t.stop() << endl;

  if (false) {//Q.size() / uniquePts.size() >= 2) {
    return concurrentHull(uniquePts);
  } else {
    return uniquePts;
  }
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc = 0) {
  using namespace std;
  using namespace parlay;
  using pt = pointVertex;

  size_t n = P.size();
  cout << "#-points = " << n << endl;

  timer t; t.start();

  sequence<pt> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = pt(P[i].coords());});

  sequence<pt> Q2 = concurrentHull(Q, numProc);

  cout << "> concurrent-hull-time = " << t.get_next() << endl;

  auto origin = pointOrigin();
  auto finalLinkedHull = new _hull<linkedFacet3d<pt>, pt, pointOrigin>(make_slice(Q2), origin, false);
  incrementHull3dSerial<linkedFacet3d<pt>, pt, pointOrigin>(finalLinkedHull); // todo parallelize

  cout << "> merge-hull-time = " << t.stop() << endl;

  auto out = sequence<facet3d<pargeo::fpoint<3>>>();
  finalLinkedHull->getHull<pargeo::fpoint<3>>(out);

  cout << "hull-size = " << out.size() << endl;

  delete finalLinkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
