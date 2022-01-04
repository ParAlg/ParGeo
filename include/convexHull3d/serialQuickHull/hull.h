#pragma once

#include "convexHull3d/vertex.h"
#include "convexHull3d/facet.h"
#include "convexHull3d/serialQuickHull/hullImpl.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {
  namespace hull3d {
    namespace serialQuickHull {

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(parlay::slice<pointT*, pointT*> P);

      template<class pointT>
      parlay::sequence<pargeo::hull3d::vertex<pargeo::hull3d::serialQuickHull::linkedFacet<pointT>, pointT>>
      computeVertex(parlay::slice<pointT*, pointT*> P);

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(hullTopology<pointT>* linkedHull);

    } // End namespace serialQuickHull
  } // End namespace hull3d
} // End namespace pargeo
