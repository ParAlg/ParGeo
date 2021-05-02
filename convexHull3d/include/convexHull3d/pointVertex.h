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
};

static std::ostream& operator<<(std::ostream& os, const pointVertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

class pointOrigin;

template <class vertexT>
struct linkedFacet3d {

  //#ifdef SERIAL
  //typedef vector<vertexT> seqT;
  //#else
  typedef sequence<vertexT> seqT;
  //#endif

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

#ifdef RESERVE
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

  void reassign(seqT *_seeList, pointOrigin* o) {
    delete seeList;
    seeList = _seeList;
  }


  size_t numVisiblePts() { return seeList->size(); }

  vertexT& visiblePts(size_t i) { return seeList->at(i); }

  size_t numPts() { return seeList->size(); }

  vertexT& pts(size_t i) { return seeList->at(i); }

  void clear() { seeList->clear(); }

  void push_back(vertexT v, pointOrigin* o) { seeList->push_back(v); }

  vertexT furthestSerial() {
    auto apex = vertexT();
    typename vertexT::floatT m = apex.attribute.numericKnob;
    for (size_t i=0; i<numVisiblePts(); ++i) {
      auto m2 = pargeo::signedVolumeX6(a, b, c, visiblePts(i));
      if (m2 > m) {
	m = m2;
	apex = visiblePts(i);
      }
    }
    return apex;
  }

  vertexT furthestParallel() {
    auto apex = vertexT();
    typename vertexT::floatT m = apex.attribute.numericKnob;
    apex = parlay::max_element(seeList->cut(0, seeList->size()),
			       [&](vertexT aa, vertexT bb) {
				 return pargeo::signedVolumeX6(a, b, c, aa) <
				   pargeo::signedVolumeX6(a, b, c, bb);
			       });
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > _a.attribute.numericKnob)//numericKnob)
      swap(b, c);

    seeList = new seqT();

#ifdef RESERVE
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

class pointOrigin {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;
  using facetT = linkedFacet3d<pointVertex>;
  using vertexT = pointVertex;

  vertexT o;

public:

  void setOrigin(vertexT _o) { o = _o; }

  inline vertexT get() {return o;};

  pointOrigin() {}

  template<class pt>
  pointOrigin(pt _o) {
    o = vertexT(_o.coords());
  }

  inline bool visible(facetT* f, vertexT p) {
    return pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob;
  }

  inline bool keep(facetT* f, vertexT p) {
    if (pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob)
      return f->a != p && f->b != p && f->c != p;
    else
      return false;
  }

};
