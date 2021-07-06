#pragma once

#include "convexHull3d/facet.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "parlay/sequence.h"

namespace pargeo {
  namespace hull3d {
    namespace sampling {

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(parlay::slice<pointT*, pointT*>, double fraction = 0.001);

      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      random(parlay::slice<pointT*, pointT*>, double);

      template<class pointT>
      parlay::sequence<facet<pointT>>
      grid(parlay::slice<pointT*, pointT*>, double);

      template<class pointT>
      parlay::sequence<facet<pointT>>
      projection(parlay::slice<pointT*, pointT*>, double);

      template<class pointT>
      double test(parlay::slice<pointT*, pointT*>, double);

    } // End namespace pargeo
  } // End namespace hull3d
} // End namespace pargeo
