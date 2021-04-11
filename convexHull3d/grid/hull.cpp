#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hull.h"
#include "grid.h"

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;

  size_t n = P.size();
  cout << "#-points = " << n << endl;

  timer t; t.start();

  // Estimate the grid size

  std::atomic<floatT> extrema[6];
  parallel_for(0, P.size(), [&](size_t i){
			      write_max(&extrema[0], P[i][0], std::less<floatT>());
			      write_min(&extrema[1], P[i][0], std::less<floatT>());
			      write_max(&extrema[2], P[i][1], std::less<floatT>());
			      write_min(&extrema[3], P[i][1], std::less<floatT>());
			      write_max(&extrema[4], P[i][2], std::less<floatT>());
			      write_min(&extrema[5], P[i][2], std::less<floatT>());
			    });
  floatT gsize = max(extrema[4]-extrema[5],
		     max((extrema[2]-extrema[3]),(extrema[1]-extrema[0]))) / 10;

  cout << "estimation-time = " << t.get_next() << endl;

  auto gridPoints = sequence<gridpt3d>(P.size());
  parallel_for(0, P.size(), [&](size_t i){
			      gridPoints[i] = gridpt3d(P[i].coords());
			      gridPoints[i].attribute = gridAtt3d(P[i], gsize);
			    });

  hashtable<hash_numeric<size_t>> table(P.size(), hash_numeric<size_t>());

  auto sample =
    parlay::filter(make_slice(gridPoints), [&](gridpt3d p) {
				     if (table.insert(p.attribute.id)) return true;
				     else return false;
				   });

  cout << "#-grids = " << sample.size() << endl;

  cout << "compute-grid-time = " << t.get_next() << endl;

  sequence<facet3d<gridpt3d>> H = incrementHull3dSerial<gridpt3d>(make_slice(sample));
  //auto H = incrementHull3d<pargeo::fpoint<3>>(make_slice(P));

  cout << "compute-hull-time = " << t.get_next() << endl;

  cout << "hull-size = " << H.size() << endl;

  // Returns true if the box intersect with hull
  auto boxCrossHull = [&](gridpt3d* box) {
			// 0: in
			// 1: out
			int ok = 0;
			for (auto f: H) {
			  if (visibleCast(&f, box[0])) {
			    ok = 1;
			    break;}
			}

			for (int i=1; i<8; ++i) {
			  gridpt3d b = box[i];
			  for (auto f: H) {
			    if (visibleCast(&f, b)) {
			      ok ^= 1;
			      break;}
			  }
			  if (ok == 1) return true;
			}

			return false;
		      };

  sequence<size_t> grids = table.entries();

  hashtable<hash_numeric<size_t>> table2(grids.size(), hash_numeric<size_t>());

  auto remainingBox =
    parlay::filter(make_slice(grids), [&](size_t g) {
					gridpt3d box[8];
					auto cs = unpack(g);
					//todo simplify
					box[0][0] = gsize * get<0>(cs);
					box[0][1] = gsize * get<1>(cs);
					box[0][2] = gsize * get<2>(cs);
					box[1] = box[0]; box[1][0]+=gsize;
					box[2] = box[0]; box[2][1]+=gsize;
					box[3] = box[0]; box[3][2]+=gsize;
					box[4] = box[0]; box[4][0]+=gsize; box[4][1]+=gsize;
					box[5] = box[0]; box[5][1]+=gsize; box[5][2]+=gsize;
					box[6] = box[0]; box[6][2]+=gsize; box[6][0]+=gsize;
					box[7] = box[0]; box[7][0]+=gsize; box[7][1]+=gsize; box[7][2]+=gsize;
					if (boxCrossHull(box)) {
					  table2.insert(g);
					  return true;
					} else {
					  return false;
					}
				      });


  auto sample2 =
    parlay::filter(make_slice(gridPoints), [&](gridpt3d p) {
					     if (table2.find(p.attribute.id) != -1) return true;
					     else return false;
					   });

  cout << "filtering-time = " << t.get_next() << endl;

  cout << "#-grids = " << table2.count() << endl;
  cout << "#-rem-points = " << sample2.size() << endl;

  sequence<facet3d<gridpt3d>> H2 = incrementHull3dSerial<gridpt3d>(make_slice(sample2));
  cout << H2.size() << endl;

  cout << "hull-time = " << t.get_next() << endl;

  return sequence<facet3d<pargeo::fpoint<3>>>(); //todo
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
