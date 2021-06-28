#pragma once

#include "facet.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
    hull3dSearch(parlay::sequence<pargeo::fpoint<3>> &);

}
