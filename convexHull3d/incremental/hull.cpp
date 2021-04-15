#include "parlay/hash_table.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hullTopology.h"
#include "pointVertex.h"
#include "hull.h"

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();
  cout << "#-points = " << n << endl;
  cout << "#-procs = " << num_workers() << endl;

  sequence<pointVertex> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = pointVertex(P[i].coords());
			    });

  // Create an initial simplex
  auto linkedHull = new _hull<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(make_slice(Q));

  incrementHull3d<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(linkedHull, numProc);

  // linkedHull is translated
  // getHull will undo the translation
  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

  cout << out.size() << endl;

  delete linkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &, size_t);
