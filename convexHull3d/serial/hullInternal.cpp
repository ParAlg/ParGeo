#include "convexHull3d/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "serialHull.h"
#include "incremental.h"
#include "vertex.h"

parlay::sequence<vertex>
pargeo::hullInternal::hull3dSerialInternal(parlay::slice<vertex*, vertex*> Q) {
  using namespace std;
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;

  auto origin = pointOrigin();

  auto linkedHull = new serialHull<linkedFacet3d<vertex>, vertex, pointOrigin>(Q, origin);

  incrementHull3dSerial<linkedFacet3d<vertex>, vertex>(linkedHull);

  return linkedHull->getHullVertices<vertex>();
}
