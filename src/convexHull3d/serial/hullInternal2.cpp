#include "convexHull3d/serialHull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "incremental.h"
#include "serialHull.h"
#include "vertex.h"

parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
pargeo::hullInternal::hull3dSerialInternal2(parlay::slice<pargeo::hullInternal::vertex*, pargeo::hullInternal::vertex*> Q) {
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;
  using vertexT = pargeo::hullInternal::vertex;

  auto origin = pointOrigin();

  auto linkedHull = new serialHull<linkedFacet3d<vertexT>, vertexT, pointOrigin>(make_slice(Q), origin);

  incrementHull3dSerial<linkedFacet3d<vertexT>, vertexT>(linkedHull);

  auto out = sequence<facetT>();

  linkedHull->getHull<pointT>(out);

  delete linkedHull;

  return out;
}
