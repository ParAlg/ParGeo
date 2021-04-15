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

  size_t l = 7;

  auto tree = gridTree3d<pointT>(make_slice(P), l+1);
  cout << "build-grid-time = " << t.get_next() << endl;

  sequence<gridVertex> Q = tree.level(l);

  cout << "extract-level-time = " << t.get_next() << endl;

  // ofstream myfile;
  // myfile.open("hull.txt", std::ofstream::trunc);
  // for (auto q: Q) {
  //   myfile << q[0] << " " << q[1] << " " << q[2] << endl;
  // }
  // myfile.close();

  cout << "level-size = " << tree.levelSize(l) << endl;
  cout << "box-size = " << tree.boxSize(l) << endl;
  cout << "max-span = " << tree.span() << endl;
  cout << "approx-factor = " << sqrt(3)*tree.boxSize(l)/tree.span() << endl;

  // Create a coarse simplex
  auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(make_slice(Q));
  linkedHull->origin.setMin(tree.getMin());
  linkedHull->origin.setGridSize(tree.boxSize(l));

  incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex, gridOrigin>(linkedHull);

  cout << "hull-time = " << t.get_next() << endl;

  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

  cout << out.size() << endl;

  delete linkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
