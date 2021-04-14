#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
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

  sequence<sequence<facet3d<pt>>> subHulls(numProc);

  parallel_for(0, numProc, [&](size_t i) {
			     size_t s = i * blkSize;
			     size_t e = min(P.size(), (i+1) * blkSize);
			     //cout << s << "--" << e << endl;
			     auto linkedHull = incrementHull3dSerial(P.cut(s, e));
			     linkedHull->getHull<pt>(subHulls[i]);
			     //cout << subHulls[i].size() << endl;
			   }, 1);

  cout << "subhull-time = " << t.get_next() << endl;

  sequence<facet3d<pt>> mergeHull = parlay::flatten(subHulls);
  //cout << mergeHull.size() << endl;

  sequence<pt> hullPoints(mergeHull.size()*3);

  parallel_for(0, mergeHull.size(), [&](size_t i) {
				      hullPoints[i*3] = mergeHull[i].a;
				      hullPoints[i*3+1] = mergeHull[i].b;
				      hullPoints[i*3+2] = mergeHull[i].c;
				    });

  parlay::sort_inplace(hullPoints);

  sequence<pt> uniquePts = parlay::unique(hullPoints);

  cout << "merge-time = " << t.get_next() << endl;

  cout << "hull-2-input-size = " << uniquePts.size() << endl;
  auto finalLinkedHull = incrementHull3dSerial(make_slice(uniquePts)); // todo parallelize

  cout << "hull2-time = " << t.stop() << endl;

  auto out = sequence<facetT>();
  finalLinkedHull->getHull<pt>(out);

  cout << out.size() << endl;

  delete finalLinkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
