#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "geometry/point.h"
#include "common/get_time.h"
#include "incremental.h"
#include "hullTopology.h"
#include "hull.h"
#include "gridTree.h"

#include <iostream>
#include <fstream>

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();

  timer t;
  t.start();

  auto tree = gridTree3d<pointT>(make_slice(P));
  cout << "build-grid-time = " << t.get_next() << endl;

  size_t L = tree.numLevels();
  size_t s = 4;
  size_t e = L+1;
  sequence<gridVertex> Q = tree.level(s);

  cout << "extract-level-time = " << t.get_next() << endl;

  // {
  //   ofstream myfile;
  //   myfile.open("point.txt", std::ofstream::trunc);
  //   for (auto q: P) {
  //     myfile << q[0] << " " << q[1] << " " << q[2] << endl;
  //   }
  //   myfile.close();
  // }

  sequence<gridVertex> pts;
  sequence<gridVertex> hullPts;

  auto out = sequence<facetT>();

  for (size_t l = s; l < e; ++ l) {
    cout << ">>> processing level " << l << endl;

    cout << "num-pts = " << Q.size() << endl;
    cout << "level-size = " << tree.levelSize(l) << endl;
    cout << "box-size = " << tree.boxSize(l) << endl;
    cout << "max-span = " << tree.span() << endl;
    cout << "approx-factor = " << sqrt(3)*tree.boxSize(l)/tree.span() << endl;

    // Create a coarse simplex
    auto origin = gridOrigin();
    origin.setMin(tree.getMin());
    origin.setBoxSize(tree.boxSize(l));
    auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(make_slice(Q), origin, false);

    incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(linkedHull);

    cout << "hull-time = " << t.get_next() << endl;

    out = sequence<facetT>(); // do we need this?
    linkedHull->getHull<pointT>(out);
    hullPts = linkedHull->getHullPts<gridVertex>(); //todo remove

    pts = linkedHull->getHullVertices(); //todo remove
    cout << "hull1-size = " << pts.size() << endl;

    auto keep = sequence<size_t>(tree.levelSize(l)+1, 0);
    parallel_for(0, pts.size(),
		 [&](size_t i) {
		   keep[pts[i].attribute.i] = 1;
		 });
    Q = tree.refine(l, keep); // includes the hull points todo memory?
    cout << "refined-pts = " << Q.size() << endl;

    cout << "hull size = " << out.size() << endl;
    cout << endl;
    delete linkedHull;
  }

  // {
  //   ofstream myfile;
  //   myfile.open("hull.txt", std::ofstream::trunc);
  //   for (auto q: hullPts) {
  //     myfile << q[0] << " " << q[1] << " " << q[2] << endl;
  //   }
  //   myfile.close();
  // }

  // {
  //   ofstream myfile;
  //   myfile.open("other.txt", std::ofstream::trunc);
  //   for (auto q: Q) {
  //     myfile << q[0] << " " << q[1] << " " << q[2] << endl;
  //   }
  //   myfile.close();
  // }
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
