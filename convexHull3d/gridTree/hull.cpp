#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hullTopology.h"
#include "gridVertex.h"
#include "hull.h"
#include "gridTree.h"

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();

  auto tree = gridTree3d<pointT>(make_slice(P), 5);

  // Now bridge the two algorithms
  sequence<gridVertex> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = gridVertex(P[i].coords());
			      // Initialize meta data related to the data type
			      // (gridVertex) todo
			    });

  // Create an initial simplex
  auto linkedHull = new _hull<linkedFacet3d<gridVertex>, gridVertex>(make_slice(Q));

  incrementHull3dSerial<gridVertex>(linkedHull);

  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

  cout << out.size() << endl;

  delete linkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
