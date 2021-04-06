#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hull.h"

#include <iostream>
#include <fstream>

using namespace std;

inline tuple<short,short,short> unpack(size_t g) {
  return make_tuple((short)g, (short)(g>>16), (short)(g>>32));
}

inline size_t pack(pargeo::fpoint<3> p, pargeo::fpoint<3>::floatT gsize) {
  size_t g = 0;
  for(int i=0; i<3; ++i) {
    unsigned short tmp = floor(p[i] / gsize);
    size_t tmp2 = 0;
    tmp2 += tmp;
    g |= tmp2 << (16*i);
  }
  return g;
}

struct gridAtt3d {
  size_t id;

  gridAtt3d() {}

  gridAtt3d(pargeo::fpoint<3> p, pargeo::fpoint<3>::floatT gsize) {
    id = pack(p, gsize);
  }
};

using gridpt3d = pargeo::_point<3, pointT::floatT, pointT::floatT, gridAtt3d>;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;

  size_t n = P.size();
  cout << "#-points = " << n << endl;

  // {
  //   ofstream myfile;
  //   myfile.open("hull.txt", std::ofstream::trunc); myfile.close();
  //   myfile.open("point.txt", std::ofstream::trunc);
  //   for(size_t p=0; p<P.size(); ++p)
  //     myfile << P[p] << p << endl;
  //   myfile.close();
  // }

  // auto HH = incrementHull3dSerial<pargeo::fpoint<3>>(make_slice(P));
  // cout << HH.size() << endl;
  // return parlay::sequence<facet3d<pargeo::fpoint<3>>>();

  timer t; t.start();

  floatT gsize = 300; // todo tune

  hashtable<hash_numeric<size_t>> table(P.size(), hash_numeric<size_t>());
  auto gridPoints = sequence<gridpt3d>(P.size());
  parallel_for(0, P.size(), [&](size_t i){
			      gridPoints[i] = gridpt3d(P[i].coords());
			      gridPoints[i].attribute = gridAtt3d(P[i], gsize);
			    });

  auto sample =
    parlay::filter(make_slice(gridPoints), [&](gridpt3d p) {
				     if (table.insert(p.attribute.id)) return true;
				     else return false;
				   });

  cout << "#-sample-1 = " << sample.size() << endl;

  cout << "sampling-1-time = " << t.get_next() << endl;

  sequence<facet3d<gridpt3d>> H = incrementHull3dSerial<gridpt3d>(make_slice(sample));
  //auto H = incrementHull3d<pargeo::fpoint<3>>(make_slice(P));

  cout << "hull-1-time = " << t.get_next() << endl;

  cout << "hull-1-size = " << H.size() << endl;

  // Put the points back, but only part of them
  auto boxCrossHull = [&](gridpt3d* box) {
			// find out using H and box
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

  // ofstream myfile2;
  // myfile2.open("other.txt", std::ofstream::trunc);

  auto remainingBox =
    parlay::filter(make_slice(grids), [&](size_t g) {
					gridpt3d box[8];
					auto cs = unpack(g);
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
					  // for(int i=0; i<8; ++i)
					  //   myfile2 << box[i][0] << " " << box[i][1] << " " << box[i][2] << endl;
					  return true;
					} else {
					  return false;
					}
				      });


  // myfile2.close();
  cout << "#-grids-remain = " << table2.count() << endl;

  auto sample2 =
    parlay::filter(make_slice(gridPoints), [&](gridpt3d p) {
					     if (table2.find(p.attribute.id) != -1) return true;
					     else return false;
					   });

 // {
 //    ofstream myfile;
 //    myfile.open("hull.txt", std::ofstream::trunc);
 //    for(auto x: sample2)
 //      myfile << x[0] << " " << x[1] << " " << x[2] << endl;
 //    myfile.close();
 //  }
  cout << "#-points-remain = " << sample2.size() << endl;

  cout << "grid-2-time = " << t.get_next() << endl;

  sequence<facet3d<gridpt3d>> H2 = incrementHull3dSerial<gridpt3d>(make_slice(sample2));
  cout << H2.size() << endl;

  cout << "hull-2-time = " << t.get_next() << endl;

  return sequence<facet3d<pargeo::fpoint<3>>>(); //todo
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
