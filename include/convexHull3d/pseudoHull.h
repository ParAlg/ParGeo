#pragma once

#include "facet.h"
#include "vertex.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dPseudo(parlay::sequence<pargeo::fpoint<3>> &);

}
