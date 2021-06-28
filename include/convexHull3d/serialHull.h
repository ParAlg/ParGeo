#pragma once

#include "facet.h"
#include "vertex.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  namespace hullInternal {

    parlay::sequence<vertex>
    hull3dSerialInternal1(parlay::slice<vertex*, vertex*>);

    parlay::sequence<facet3d<fpoint<3>>>
    hull3dSerialInternal2(parlay::slice<vertex*, vertex*>);

    parlay::sequence<fpoint<3>>
    hull3dSerialInternal3(parlay::sequence<fpoint<3>> &);

  } // End namespace hullInternal

  parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
  hull3dSerial(parlay::sequence<pargeo::fpoint<3>> &);

} // End namespace pargeo
