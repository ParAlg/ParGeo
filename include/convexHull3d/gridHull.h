#pragma once

#include "facet.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
  hull3dGrid(parlay::sequence<pargeo::fpoint<3>> &, size_t s = 4, bool write = false);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dGridConcurrent(parlay::sequence<pargeo::fpoint<3>> &P, size_t s, size_t numProc);

} // End namespace pargeo
