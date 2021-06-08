#pragma once

#include "pargeo/point.h"
#include "pargeo/algebra.h"
#include "parlay/primitives.h"

// todo namespace

class vertexAtt;

using vertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt>;

template <class vertexT> struct linkedFacet3d;

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

namespace pargeo {

  namespace hullInternal {

    // Clockwise oriented 3d facet
    template <class pt, class att>
    struct _facet3d {
      using pointT = pt;
      pt a, b, c;// Indices into P
      att attribute;
      _facet3d(pt _a, pt _b, pt _c): a(_a), b(_b), c(_c) {
	if (pargeo::determinant3by3(a, b, c) > 0)
	  std::swap(b, c);
      }
    };

    parlay::sequence<vertex>
    hull3dSerialInternal(parlay::slice<vertex*, vertex*>);

  }

  template <class pt>
  using facet3d = hullInternal::_facet3d<pt, pargeo::_empty>;

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dSerial(parlay::sequence<pargeo::fpoint<3>> &);

  parlay::sequence<pargeo::fpoint<3>>
  hull3dSerialInternal(parlay::sequence<pargeo::fpoint<3>> &);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dSerialInternal(parlay::slice<vertex*, vertex*>);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dIncremental(parlay::sequence<pargeo::fpoint<3>> &, size_t numProc = 0);

  parlay::sequence<facet3d<fpoint<3>>>
  hull3dIncrementalInternal(parlay::slice<vertex*, vertex*>, size_t numProc = 0);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dConcurrent(parlay::sequence<pargeo::fpoint<3>> &, size_t numProc = 0);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dGrid(parlay::sequence<pargeo::fpoint<3>> &, size_t, bool);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dGridConcurrent(parlay::sequence<pargeo::fpoint<3>> &, size_t s = 4, size_t numProc = 0);

  parlay::sequence<facet3d<pargeo::fpoint<3>>>
  hull3dPseudo(parlay::sequence<pargeo::fpoint<3>> &);

}; // End namespace
