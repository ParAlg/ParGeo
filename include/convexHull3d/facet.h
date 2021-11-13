#pragma once

#include "pargeo/algebra.h"

namespace pargeo {
  namespace hull3d {

    /* Clockwise oriented 3d facet */

    template <typename _pointT> class facet {
    public:
      using pointT = _pointT;

      pointT a, b, c;

      facet(pointT _a, pointT _b, pointT _c): a(_a), b(_b), c(_c) {

	if (pargeo::determinant3by3(a, b, c) > 0)
	  std::swap(b, c);

      }

    };

  } // End namespace hull3d
} // End namespace pargeo
