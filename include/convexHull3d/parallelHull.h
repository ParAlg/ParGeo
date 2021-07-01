#pragma once

#include "facet.h"
#include "vertex.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
  hull3dParallel(parlay::sequence<pargeo::fpoint<3>> &, size_t numProc = 0);

}; // End namespace pargeo
