#include "convexHull3d/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "incremental.h"
#include "serialHull.h"
#include "vertex.h"

parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
pargeo::hull3dSerialInternal(parlay::slice<vertex*, vertex*> Q) {
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;

  auto origin = pointOrigin();

  auto linkedHull = new serialHull<linkedFacet3d<vertex>, vertex, pointOrigin>(make_slice(Q), origin);

  incrementHull3dSerial<linkedFacet3d<vertex>, vertex>(linkedHull);

  auto out = sequence<facetT>();

  linkedHull->getHull<pointT>(out);

  delete linkedHull;

  return out;
}
