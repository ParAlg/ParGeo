#pragma once

#include "incremental.h"

class vertexAtt;

using gridVertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt>;

class vertexAtt {
public:
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;
  size_t i;
  linkedFacet3d<gridVertex> *seeFacet;

  vertexAtt() {}

  /* Signed volume (x6) of an oriented tetrahedron (example below is positive).
       d
       o
      /|\
     / | o b
    o--o/
     a   c

     x
     origin
  */
  template <class pt>
  inline typename pt::floatT signedVolume(pt a, pt b, pt c, pt d) {
    return (a-d).dot(crossProduct3d(b-a, c-a));
  }

  // area is precomputed from oriented triang a,b,c
  template <class pt>
  inline typename pt::floatT signedVolume(pt a, pt d, pt area) {
    return (a-d).dot(area);
  }

  template<class facetT, class vertexT>
  bool visible(facetT* f, vertexT p) {
    // if (signedVolume(f->a, f->b, f->c, p) > numericKnob)
    if (signedVolume(f->a, p, f->area) > numericKnob)
      return true;
    else
      return false;
  }

  template<class facetT, class vertexT>
  bool visibleNoDup(facetT* f, vertexT p) {
    if (signedVolume(f->a, p, f->area) > numericKnob)
      return true && f->a != p && f->b != p && f->c != p;
    else
      return false;
  }

};

static std::ostream& operator<<(std::ostream& os, const gridVertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}
