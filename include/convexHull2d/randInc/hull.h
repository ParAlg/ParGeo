#pragma once

#include "parlay/sequence.h"

namespace pargeo {
  namespace hull2d {
    namespace randInc {

      template<class pointT>
      parlay::sequence<size_t>
      compute(parlay::slice<pointT*, pointT*>, size_t numProc = 0);

    }
  }
}
