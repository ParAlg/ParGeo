#pragma once

#include "incremental.h"

class vertexAtt;

using pointVertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, vertexAtt>;

template <class vertexT> struct linkedFacet3d;

class vertexAtt {
public:
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

#ifdef VERBOSE
  size_t i;
#endif
  linkedFacet3d<pointVertex> *seeFacet;
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

static std::ostream& operator<<(std::ostream& os, const pointVertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template <class vertexT>
struct linkedFacet3d {
  //using vertexT = vertex3d;
#ifdef SERIAL
  typedef vector<vertexT> seqT;
#else
  typedef sequence<vertexT> seqT;
#endif

  vertexT a, b, c;
  linkedFacet3d *abFacet;
  linkedFacet3d *bcFacet;
  linkedFacet3d *caFacet;
  vertexT area;

  bool isAdjacent(linkedFacet3d* f2) {
    vertexT V[3]; V[0] = a; V[1] = b; V[2] = c;
    for (int i=0; i<3; ++i) {
      auto v1 = V[i];
      auto v2 = V[(i+1)%3];
      if ( (f2->a == v1 && f2->b == v2) || (f2->a == v2 && f2->b == v1) ) return true;
      if ( (f2->b == v1 && f2->c == v2) || (f2->b == v2 && f2->c == v1) ) return true;
      if ( (f2->c == v1 && f2->a == v2) || (f2->c == v2 && f2->a == v1) ) return true;
    }
    return false;
  }

#ifndef SERIAL
  // Stores the minimum memory address of the seeFacet of the reserving vertices
  std::atomic<size_t> reservation;

  void reserve(vertexT& p) {
    parlay::write_min(&reservation, (size_t)p.attribute.seeFacet, std::less<size_t>());
  }

  bool reserved(vertexT& p) {
    return reservation == (size_t)p.attribute.seeFacet;
  }
#endif

  seqT *seeList;

  void reassign(seqT *_seeList) {
    delete seeList;
    seeList = _seeList;
  }

  vertexT& at(size_t i) { return seeList->at(i); }

  size_t size() { return seeList->size(); }

  void clear() { seeList->clear(); }

  void push_back(vertexT v) { seeList->push_back(v); }

  vertexT furthest() {
    auto apex = vertexT();
    typename vertexT::floatT m = apex.attribute.numericKnob;

#ifdef SERIAL
    for (size_t i=0; i<size(); ++i) {
      auto m2 = at(i).attribute.signedVolume(a, b, c, at(i));
      if (m2 > m) {
	m = m2;
	apex = at(i);
      }
    }
#else
    apex = parlay::max_element(seeList->cut(0, seeList->size()),
			       [&](vertexT aa, vertexT bb) {
				 return aa.attribute.signedVolume(a, b, c, aa) <
				   bb.attribute.signedVolume(a, b, c, bb);
			       });
#endif
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > _a.attribute.numericKnob)//numericKnob)
      swap(b, c);

    seeList = new seqT();

#ifndef SERIAL
    reservation = -1; // (unsigned)
#endif

    area = crossProduct3d(b-a, c-a);
  }

  ~linkedFacet3d() {
    delete seeList;
  }

};

template<class vertexT>
static std::ostream& operator<<(std::ostream& os, const linkedFacet3d<vertexT>& v) {
#ifdef VERBOSE
  os << "(" << v.a.attribute.i << "," << v.b.attribute.i << "," << v.c.attribute.i << ")";
#else
  os << "(" << v.a << "," << v.b << "," << v.c << ")";
#endif
  return os;
}
