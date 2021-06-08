#include "convexHull3d/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "serialHull.h"
#include "incremental.h"
#include "vertex.h"

parlay::sequence<pargeo::fpoint<3>>
pargeo::hull3dSerialInternal(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace std;
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;

  sequence<vertex> Q(P.size());

  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = vertex(P[i].coords());
			    });

  auto origin = pointOrigin();

  auto linkedHull = new serialHull<linkedFacet3d<vertex>, vertex, pointOrigin>(make_slice(Q), origin);

  incrementHull3dSerial<linkedFacet3d<vertex>, vertex>(linkedHull);

  return linkedHull->getHullVertices<pointT>();
}
