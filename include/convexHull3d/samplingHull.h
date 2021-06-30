#pragma once

#include "facet.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "parlay/sequence.h"

namespace pargeo {

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dSampling(parlay::sequence<pargeo::fpoint<3>> &, double fraction = 0.001);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dRandomSampling(parlay::sequence<pargeo::fpoint<3>> &P, double fraction);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dGridSampling(parlay::sequence<pargeo::fpoint<3>> &P, double fraction);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dRandomProjection(parlay::sequence<pargeo::fpoint<3>> &P, double fraction);

  double
  testHull(parlay::sequence<pargeo::fpoint<3>> &, double);

} // End namespace pargeo
