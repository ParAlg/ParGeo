#pragma once

#include "pargeo/point.h"

template <class vertexT> struct linkedFacet3d;

namespace pargeo {
  namespace hullInternal {

    // Internal vertex

    class vertexAtt;

    using vertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt>;

    // should probably template the facet
    class vertexAtt {
    public:
      static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

      // #ifdef VERBOSE
      //   size_t i;
      // #endif
      linkedFacet3d<vertex> *seeFacet;
      vertexAtt() {}
    };

    static std::ostream& operator<<(std::ostream& os, const vertex& v) {
      for (int i=0; i<v.dim; ++i)
	os << v.x[i] << " ";
      return os;
    }

  } // End namespace hullInternal
} // End namespace pargeo
