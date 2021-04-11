#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hull.h"
#include "grid.h"

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc = 0) {
  using namespace std;
  using namespace parlay;
  using pt = pargeo::fpoint<3>;

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
			     subHulls[i] = incrementHull3dSerial(P.cut(s, e));
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
  auto H = incrementHull3dSerial(make_slice(uniquePts)); // todo parallelize

  cout << "hull2-time = " << t.stop() << endl;
  cout << H.size() << endl;

  return H;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
