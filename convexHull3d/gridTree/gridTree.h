#pragma

#include <bitset>
#include "parlay/sequence.h"
#include "geometry/point.h"

using namespace parlay;

class gridAtt3d;

using gridVertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, gridAtt3d>;

template <class vertexT> struct linkedFacet3d;

struct gridAtt3d {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;
  static constexpr size_t maxRange = 65535;
  using pt = pargeo::fpoint<3>;
  using floatT = pt::floatT;

  inline size_t shift(size_t x) {
    if (x > maxRange) {
      cout << "Grid number " << x << endl;
      throw std::runtime_error("Grid id exceed 16 bits, make min-grid size larger.");
    }

    x = (x | (x << 16)) & 0x001f0000ff0000ff;
    x = (x | (x <<  8)) & 0x100f00f00f00f00f;
    x = (x | (x <<  4)) & 0x10c30c30c30c30c3;
    x = (x | (x <<  2)) & 0x1249249249249249;
    return x;
  }

  void computeId(pargeo::fpoint<3> p, pargeo::fpoint<3> pMin, floatT gSize) {
    size_t x = floor( (p[0] - pMin[0]) / gSize );
    size_t y = floor( (p[1] - pMin[1]) / gSize );
    size_t z = floor( (p[2] - pMin[2]) / gSize );

    x = shift(x);
    y = shift(y);
    z = shift(z);
    id = x | (y << 1) | (z << 2);
  }

  // Level 0 is the coarsest level
  inline size_t getLevel(size_t l) {
    // the highest 16 bits are empty
    // starting taking 3-bit numbers
    return id >> (48 - (l+1)*3);
  }

  gridAtt3d() {}

  size_t id; // grid id

  size_t i; // index, todo check if needed

  linkedFacet3d<gridVertex> *seeFacet;
};

static std::ostream& operator<<(std::ostream& os, const gridVertex& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template<class pt>
class gridTree3d {
  using floatT = gridVertex::floatT;
  static constexpr size_t maxRange = 65535;
  static constexpr size_t maxBit = 16;

  sequence<gridVertex> P;

  sequence<size_t>** pointers;

  size_t L;

  floatT gSize;

  pt pMin;

  floatT maxSpan;
public:

  sequence<size_t>* levelIdx(size_t l) {
    return pointers[l];
  }

  floatT boxSize(size_t l) {
    return gSize * pow(2, maxBit - l);
  }

  floatT span() {
    return maxSpan;
  }

  pt getMin() {
    return pMin;
  }

  floatT levelSize(size_t l) {
    return pointers[l]->size() - 1;
  }

  sequence<gridVertex> level(size_t l) {
    sequence<gridVertex> out(pointers[l]->size()-1);
    parallel_for(0, pointers[l]->size()-1,
		 [&](size_t i){
			      out[i] = P[pointers[l]->at(i)];
			    });
    return out;
  }

  gridVertex at(size_t l, size_t i) {
    return P[pointers[l]->at(i)];
  }

  gridVertex at(size_t i) {
    return P[i];
  }

  ~gridTree3d() {
    //todo
  }

  gridTree3d(slice<pt*, pt*> _P, size_t _L) {
    L = _L;

    std::atomic<floatT> extrema[6];
    for (int i=0; i<3; ++i) {
      extrema[i*2] = _P[0][i];
      extrema[i*2+1] = _P[0][i];
    }
    parallel_for(0, _P.size(), [&](size_t i){
				write_max(&extrema[0], _P[i][0], std::less<floatT>());
				write_min(&extrema[1], _P[i][0], std::less<floatT>());
				write_max(&extrema[2], _P[i][1], std::less<floatT>());
				write_min(&extrema[3], _P[i][1], std::less<floatT>());
				write_max(&extrema[4], _P[i][2], std::less<floatT>());
				write_min(&extrema[5], _P[i][2], std::less<floatT>());
			      });
    maxSpan = max(extrema[4]-extrema[5],
		  max((extrema[2]-extrema[3]),(extrema[1]-extrema[0])));
    size_t _maxRange = maxRange;
    gSize = 1.01 * maxSpan / maxRange;

    pMin[0] = extrema[1];
    pMin[1] = extrema[3];
    pMin[2] = extrema[5];
    cout << "untranslated-pmin = " << pMin << endl;

    cout << "grid-size = " << gSize << endl;

    P = sequence<gridVertex>(_P.size());
    parallel_for(0, P.size(), [&](size_t i){
				P[i] = gridVertex(_P[i].coords());
				P[i].attribute = gridAtt3d();
				P[i].attribute.computeId(_P[i], pMin, gSize);
			      });

    parlay::sort_inplace(make_slice(P), [&](gridVertex const& a, gridVertex const& b){
					  return a.attribute.id < b.attribute.id;});

    sequence<size_t> flag(P.size()+1);

    // The the coarser L levels
    pointers = (sequence<size_t>**) malloc(sizeof(sequence<size_t>*)*L);
    for (int l=0; l<L; ++l) {
      flag[0] = 1;
      parallel_for(1, P.size(), [&](size_t i){
				  if (P[i].attribute.getLevel(l) != P[i-1].attribute.getLevel(l)) {
				    flag[i] = 1;
				  } else flag[i] = 0;
				});
      size_t numGrids = parlay::scan_inplace(flag.cut(0, P.size()));
      flag[P.size()] = numGrids;
      pointers[l] = new sequence<size_t>(numGrids+1);
      parallel_for(0, P.size(), [&](size_t i){
				  if (flag[i] != flag[i+1]) {
				    pointers[l]->at(flag[i]) = i;
				  }
				});
      pointers[l]->at(numGrids) = P.size();
      cout << "lvl " << l << " = " << numGrids << endl;
    }

  }
};

template <class vertexT>
struct linkedFacet3d {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

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
      auto m2 = pargeo::signedVolumeX6(a, b, c, at(i));
      if (m2 > m) {
	m = m2;
	apex = at(i);
      }
    }
#else
    apex = parlay::max_element(seeList->cut(0, seeList->size()),
			       [&](vertexT aa, vertexT bb) {
				 return pargeo::signedVolumeX6(a, b, c, aa) <
				   pargeo::signedVolumeX6(a, b, c, bb);
			       });
#endif
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > _a.attribute.numericKnob)
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

class gridOrigin {
  using facetT = linkedFacet3d<gridVertex>;
  using vertexT = gridVertex;
  using floatT = typename pargeo::fpoint<3>::floatT;
  static constexpr floatT numericKnob = 1e-5;

  vertexT o; // translation origin

  vertexT pMin; // translated pmin of the grid

  floatT gridSize; // gridSize of the level being processed

public:
  void setGridSize(floatT _gridSize) {
    gridSize = _gridSize;
  }

  template<class pt>
  void setMin(pt _pMin) {
    pMin = vertexT(_pMin.coords());
    pMin = pMin - o; // translate the pmin by origin (same as other points in the hull)
  }

  inline vertexT get() {return o;};

  gridOrigin() {}

  template<class pt>
  gridOrigin(pt _o) {
    o = vertexT(_o.coords());
  }

  inline bool visible(facetT* f, vertexT p) {
    if (pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob)
      return true;
    else
      return false;
  }

  inline bool visibleNoDup(facetT* f, vertexT p) {
    if (pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob)
      return true && f->a != p && f->b != p && f->c != p;
    else
      return false;
  }

};
