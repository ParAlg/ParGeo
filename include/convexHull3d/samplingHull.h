#pragma once

#include "facet.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dSampling(parlay::sequence<pargeo::fpoint<3>> &, double);

  double
  testHull(parlay::sequence<pargeo::fpoint<3>> &, double);

} // End namespace pargeo
