#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "convexHull3d/incremental.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/hull.h"
#include "gridTree.h"

#include <iostream>
#include <fstream>

//#define WRITE

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();

  pargeo::timer t;
  t.start();

  auto tree = gridTree3d<pointT>(make_slice(P));
  cout << "build-grid-time = " << t.get_next() << endl;

  size_t L = tree.numLevels();
  size_t s = 3;
  size_t e = L;
  cout << "L = " << L << endl;

  // for (size_t l = s; l < e; ++ l) {
  //   sequence<gridVertex> Q = tree.level(l);
  //   cout << "---" << endl;
  //   cout << "level-" << l << "-size = " << Q.size() << endl;
  //   cout << "level-" << l << "-box-size = " << tree.boxSize(l) << endl;
  //   cout << "extract-level-time = " << t.get_next() << endl;
  // }
  // cout << endl;

  sequence<gridVertex> Q = tree.level(s);

  auto out = sequence<facetT>();

  for (size_t l = s; l < e; ++ l) {
    bool lastRound = tree.levelSize(l) == P.size() || l == e-1;

    cout << "--- level " << l << endl;
    cout << "level-size = " << tree.levelSize(l) << endl;
    cout << "box-size = " << tree.boxSize(l) << endl;
    cout << "num-pts = " << Q.size() << endl;

#ifdef WRITE
    {
      ofstream myfile;
      myfile.open("point.txt", std::ofstream::trunc);
      for (auto q: Q) {
	myfile << q[0] << " " << q[1] << " " << q[2] << endl;
	//myfile << q[0] << " " << q[1] << " " << q[2] << " " << q.attribute.i << endl;
      }
      myfile.close();
    }
#endif

    // Create a coarse simplex
    auto origin = gridOrigin();
    origin.setMin(tree.getMin());
    if (lastRound) {
      origin.setBoxSize(-1);
    } else {
      origin.setBoxSize(tree.boxSize(l));
    }
    auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(make_slice(Q), origin, false);
    origin = linkedHull->getOrigin();

    incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(linkedHull);

    cout << "hull-time = " << t.get_next() << endl;

    auto pts = linkedHull->getHullPts();

#ifdef WRITE
    linkedHull->writeHull("facet.txt");
#endif

    // for (auto x: Q) {
    //   size_t id = x.attribute.i;
    //   size_t nextLevelId = tree.pointers[l]->at(x.attribute.i);
    //   //cout << nextLevelId << " ";
    //   cout << id << " " << nextLevelId << ": ";
    //   auto children = tree.children(l, x.attribute.i);
    //   for (auto y: children) {
    //   	cout << y.attribute.i << " ";
    //   }
    //   cout << endl;
    // }
    // cout << endl;

#ifdef WRITE
    // { // points in the next level
    //   ofstream myfile;
    //   myfile.open("blue.txt", std::ofstream::trunc);
    //   // origin.writeCorners(Q[0], myfile);
    //   // origin.writeCorners(Q[2], myfile);
    //   for (auto q: Q) {
    // 	origin.writeCorners(q, myfile);
    //   }
    //   myfile.close();
    // }
#endif

    auto keep = sequence<size_t>(tree.levelSize(l)+1, 0);
    parallel_for(0, pts.size(),
		 [&](size_t i) {
		   keep[pts[i].attribute.i] = 1;
		 });
    Q = tree.nextLevel(l, keep);
    cout << "refined-pts = " << Q.size() << endl;

#ifdef WRITE
    auto P2 = tree.level(l+1);
    { // points in the next level but not considered next round
      ofstream myfile;
      myfile.open("red.txt", std::ofstream::trunc);
      for (auto p: P2) {
        bool ok = true;
        for (auto q: Q) {
	  if (p == q) {
	    ok = false;
	    break;
	  }
        }
        if (ok) myfile << p << endl;
      }
      myfile.close();
    }
#endif

    if (lastRound) {
      linkedHull->getHull<pointT>(out);
      delete linkedHull;
      break;
    }
    delete linkedHull;
  }

  cout << "hull-size = " << out.size() << endl;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
