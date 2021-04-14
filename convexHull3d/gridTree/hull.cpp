#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
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

  auto tree = gridTree3d<pointT>(make_slice(P), 10);

  sequence<gridVertex> Q = tree.level(4);

  // ofstream myfile;
  // myfile.open("hull.txt", std::ofstream::trunc);
  // for (auto q: Q) {
  //   myfile << q[0] << " " << q[1] << " " << q[2] << endl;
  // }
  // myfile.close();

  cout << Q.size() << endl;

  // Create a coarse simplex
  auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex>(make_slice(Q));

  incrementHull3dSerial<linkedFacet3d<gridVertex>, gridVertex>(linkedHull);

  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

  cout << out.size() << endl;

  delete linkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
