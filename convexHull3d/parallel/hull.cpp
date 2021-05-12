#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "convexHull3d/incremental.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/pointVertex.h"
#include "convexHull3d/hull.h"

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3dParallel(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;

  size_t n = P.size();
#ifndef SILENT
  cout << "#-points = " << n << endl;
#endif
  sequence<pointVertex> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = pointVertex(P[i].coords());
			      // Initialize meta data related to the data type
			      // (pointVertex)
			      // Nothing here
			    });

  // Create an initial simplex
  auto origin = pointOrigin();
  auto linkedHull = new _hull<linkedFacet3d<pointVertex>, pointVertex, pointOrigin>(make_slice(Q), origin, false);

  incrementHull3dParallel<linkedFacet3d<pointVertex>, pointVertex>(linkedHull);

  // linkedHull is translated
  // getHull will undo the translation
  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

#ifndef SILENT
  cout << out.size() << endl;
#endif
  delete linkedHull;
  return out;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3dSerial(parlay::sequence<pargeo::fpoint<3>> &);
