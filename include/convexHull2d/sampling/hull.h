#pragma once

#include "parlay/sequence.h"

namespace pargeo {
  namespace hull2d {
    namespace sampling {

      template<class pointT>
      parlay::sequence<size_t>
      compute(parlay::slice<pointT*, pointT*>, double fraction = 0.001);

    }
  }
}
