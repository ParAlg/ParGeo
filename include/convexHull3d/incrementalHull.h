#pragma once

#include "facet.h"
#include "vertex.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
  hull3dIncremental(parlay::sequence<pargeo::fpoint<3>> &, size_t numProc = 0);

  parlay::sequence<facet3d<fpoint<3>>>
  hull3dIncrementalInternal(parlay::slice<pargeo::hullInternal::vertex*, pargeo::hullInternal::vertex*>, size_t numProc = 0);

}; // End namespace pargeo
