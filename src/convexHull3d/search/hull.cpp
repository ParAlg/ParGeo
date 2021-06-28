#include "convexHull3d/bruteforceHull.h"
#include "convexHull3d/searchHull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "incremental.h"
#include "serialHull.h"
#include "vertex.h"

parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
pargeo::hull3dSearch(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;

  if (P.size() < 10) {
    return hull3dBruteforce(P);
  }

  sequence<vertex> Q(P.size());

  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = vertex(P[i].coords());
			    });

  auto origin = pointOrigin();

  auto linkedHull = new serialHull<linkedFacet3d<vertex>, vertex, pointOrigin>(make_slice(Q), origin);

  incrementHull3dSerial<linkedFacet3d<vertex>, vertex>(linkedHull);

  auto out = sequence<facetT>();

  linkedHull->getHull<pointT>(out);

  delete linkedHull;

  return out;
}
