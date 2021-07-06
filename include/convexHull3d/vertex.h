#pragma once

#include "pargeo/point.h"

// template <class vertexT> struct linkedFacet3d;

namespace pargeo {
  namespace hull3d {

    template<class linkedFacet3d, class T>
    class vertex: public T {

    public:

      using pointT = T;
      using floatT = typename T::floatT;

      linkedFacet3d *seeFacet;

      vertex(floatT* p): T(p) {};

      vertex(floatT* p, typename T::attT a): T(p, a) {};

      vertex(): T() {};

      vertex<linkedFacet3d, T> operator/(floatT dv) {
	floatT xx[T::dim];
	for (int i=0; i<T::dim; ++i) xx[i] = T::x[i]/dv;
	return vertex<linkedFacet3d, T>(xx, T::attribute);}
    };

    // Internal vertex

    // template<class linkedFacet3d>
    // class vertexAtt;

    // template<class linkedFacet3d>
    // using vertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt<linkedFacet3d>>;

    // template<class linkedFacet3d>
    // class vertexAtt {
    // public:
    //   static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

    //   // #ifdef VERBOSE
    //   //   size_t i;
    //   // #endif
    //   linkedFacet3d *seeFacet;
    //   vertexAtt() {}
    // };

    // template<class FT, class PT>
    // static std::ostream& operator<<(std::ostream& os, const vertex<FT, PT>& v) {
    //   for (int i=0; i<v.dim; ++i)
    // 	os << v.x[i] << " ";
    //   return os;
    // }

  } // End namespace hullInternal
} // End namespace pargeo
