#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "pointVertex.h"
#include "hullTopology.h"
#include "incremental.h"
#include "hull.h"

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc = 0) {
  using namespace std;
  using namespace parlay;
  using pt = pargeo::fpoint<3>;
  using floatT = pargeo::fpoint<3>::floatT;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();
  cout << "#-points = " << n << endl;

  timer t; t.start();

  if (!numProc) numProc = num_workers();

  numProc *= 8;

  size_t blkSize = floor(P.size() / numProc);

  while (blkSize < 100) {
    numProc -= 1;
    blkSize = floor(P.size() / numProc);
  }

  cout << "#subproblems = " << numProc << endl;

  sequence<sequence<pt>> subHulls(numProc);

  sequence<pointVertex> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = pointVertex(P[i].coords());});

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = min(P.size(), (i+1) * blkSize);
			     //cout << s << "--" << e << endl;
			     auto origin = pointOrigin();
			     auto linkedHull = new _hull<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(Q.cut(s, e), origin);
			     incrementHull3dSerial<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(linkedHull);
			     subHulls[i] = linkedHull->getHullPts<pt>();
			   }, 1);

  cout << "subhull-time = " << t.get_next() << endl;

  sequence<pt> uniquePts = parlay::flatten(subHulls);

  cout << "merge-time = " << t.get_next() << endl;

  cout << "hull-2-input-size = " << uniquePts.size() << endl;

  sequence<pointVertex> Q2(uniquePts.size());

  parallel_for(0, uniquePts.size(), [&](size_t i) {
			      Q2[i] = pointVertex(uniquePts[i].coords());});

  auto origin = pointOrigin();
  auto finalLinkedHull = new _hull<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(make_slice(Q2), origin);
  incrementHull3dSerial<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(finalLinkedHull); // todo parallelize

  cout << "hull2-time = " << t.stop() << endl;

  auto out = sequence<facetT>();
  finalLinkedHull->getHull<pt>(out);

  cout << out.size() << endl;

  delete finalLinkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
