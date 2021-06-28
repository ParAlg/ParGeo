#pragma once

#include "pargeo/algebra.h"

namespace pargeo {

  // Clockwise oriented 3d facet
  template <class pt>
  struct facet3d {
    using pointT = pt;
    pt a, b, c;// Indices into P
    facet3d(pt _a, pt _b, pt _c): a(_a), b(_b), c(_c) {
      if (pargeo::determinant3by3(a, b, c) > 0)
	std::swap(b, c);
    }
  };

} // End namespace pargeo
