#pragma once

#include "convexHull3d/facet.h"
// #include "convexHull3d/serialQuickHull/hullImpl.h"
#include "convexHull3d/vertex.h"
#include "pargeo/point.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"

namespace pargeo {
  namespace hull3d {
    namespace serialQuickHull {

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(parlay::slice<pointT*, pointT*> P);

    } // End namespace serialQuickHull
  } // End namespace hull3d
} // End namespace pargeo
