#pragma once

namespace pargeo {
  namespace hull3d {
    namespace divideConquer {
      template<class pointT>
      parlay::sequence<pargeo::hull3d::facet<pointT>>
      compute(parlay::slice<pointT*, pointT*>, size_t numProc = 0);
    }
  }
}
