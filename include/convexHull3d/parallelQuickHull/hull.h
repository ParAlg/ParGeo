#pragma once

#include "convexHull3d/facet.h"
#include "convexHull3d/vertex.h"
#include "convexHull3d/parallelQuickHull/hullImpl.h"

#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

namespace pargeo {
  namespace hull3d {
    namespace parallelQuickHull {

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT> >
      compute(parlay::slice<pointT*, pointT*>, size_t = 0, bool randomized = true);

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(parlay::slice<
	      pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pointT>, pointT>*,
	      pargeo::hull3d::vertex<pargeo::hull3d::parallelQuickHull::linkedFacet<pointT>, pointT>*
	      >, size_t = 0);

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(hullTopology<pointT>*, size_t = 0);

    } // End namespace parallelQuickHull
  } // End namespace hull3d
}; // End namespace pargeo
