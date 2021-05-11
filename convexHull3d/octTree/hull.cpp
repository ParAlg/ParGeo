#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "convexHull3d/incremental.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/hull.h"
#include "gridVertex.h"
#include "octTree.h"
#include <iostream>
#include <fstream>

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3dGrid(parlay::sequence<pargeo::fpoint<3>> &P, size_t s = 0, size_t skip = 1, bool write = false) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();

  pargeo::timer t;
  pargeo::timer tr;
  t.start();

  // auto tree = octTree<gridVertex>(make_slice(P));

  auto levels = sequence<size_t>();
  //levels.push_back(15);
  levels.push_back(3);
  //levels.push_back(2);
  auto tree = octTree<gridVertex>(make_slice(P), levels);

  cout << ">>> build-grid-time = " << t.get_next() << endl;

  size_t L = tree.numLevels();

  sequence<gridVertex> Q = tree.level(0); // Get the coarsest level points
  sequence<size_t> keep;

  if (write) {
    ofstream myfile;
    myfile.open("point.txt", std::ofstream::trunc);
    for (auto q: P) {
      myfile << q[0] << " " << q[1] << " " << q[2] << endl;
      //myfile << q[0] << " " << q[1] << " " << q[2] << " " << q.attribute.i << endl;
    }
    myfile.close();
  }

  for (size_t l = 0; l < tree.numLevels(); ++l) {
    tr; tr.start();

    cout << "--------------- level " << l << endl;

    if (l > 0)
      Q = tree.nextLevel(l-1, keep);
    cout << " " << Q.size() << " / " << tree.levelSize(l) << endl;

    if (write) {
      ofstream myfile;
      myfile.open("blue.txt", std::ofstream::trunc);
      for (auto q: Q) {
	myfile << q[0] << " " << q[1] << " " << q[2] << endl;
	//myfile << q[0] << " " << q[1] << " " << q[2] << " " << q.attribute.i << endl;
      }
      myfile.close();
    }

    // Create a coarse simplex

    auto origin = gridOrigin();

    origin.setMin(tree.getMin());

    origin.setBoxSize(tree.boxSize(l));

    auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(make_slice(Q), origin, true);
    origin = linkedHull->getOrigin();

    incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(linkedHull);

    auto pts = linkedHull->getHullPts();

    keep = sequence<size_t>(tree.levelSize(l) + 1, 0);
    parallel_for(0, pts.size(),
		 [&](size_t i) {
		   keep[pts[i].attribute.i] = 1;
		 });

    if (write) {
      linkedHull->writeHull("facet.txt");

      auto P2 = tree.level(l+1);
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

    delete linkedHull;
    cout << " round-time = " << tr.stop() << endl;

  }

  auto finalQ = tree.getPoints<pointT>(tree.numLevels()-1, keep);
  cout << "--------------- final level " << endl;
  tr.start();
  cout << " " << finalQ.size() << " / " << P.size() << endl;
  auto result = hull3d(finalQ);
  cout << " round-time = " << tr.stop() << endl;

  cout << ">>> total-hull-time = " << t.get_next() << endl;

  cout << "hull size = " << result.size() << endl;
  return result;

  // cout << "hull-size = " << out.size() << endl;
  // return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &, size_t, size_t);
