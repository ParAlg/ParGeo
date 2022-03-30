#pragma once

#include "parlay/sequence.h"

namespace pargeo {

  template<typename pt>
  parlay::sequence<pt> zorderSort2d(parlay::sequence<pt>& P);

  template<typename pt>
  parlay::sequence<pt> zorderSort3d(parlay::sequence<pt>& P);

  template<typename pt>
  void zorderSortInPlace2d(parlay::sequence<pt>& P);

  template<typename pt>
  void zorderSortInPlace3d(parlay::sequence<pt>& P);

} // End namespace

#include "mortonSortImpl.h"
