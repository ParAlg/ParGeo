#pragma once

#include "pargeo/algebra.h"

namespace pargeo {
  namespace hull3d {

    // Clockwise oriented 3d facet
    template <class pt>
    struct facet {
      using pointT = pt;
      pt a, b, c;// Indices into P
      facet(pt _a, pt _b, pt _c): a(_a), b(_b), c(_c) {
	if (pargeo::determinant3by3(a, b, c) > 0)
	  std::swap(b, c);
      }
    };

  } // End namespace hull3d
} // End namespace pargeo
