#pragma

#include <iostream>
#include <fstream>
#include <tuple>
#include <bitset>
#include "parlay/sequence.h"
#include "pargeo/point.h"

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

  void computeId(pargeo::fpoint<3> p, pargeo::fpoint<3> pMin, floatT boxSize) {
    size_t x = floor( (p[0] - pMin[0]) / boxSize );
    size_t y = floor( (p[1] - pMin[1]) / boxSize );
    size_t z = floor( (p[2] - pMin[2]) / boxSize );

    x = shift(x);
    y = shift(y);
    z = shift(z);
    id = x | (y << 1) | (z << 2);
  }

  inline size_t gridId(size_t l) {
    // The highest 16 bits are empty, each level is 3-bits
    // <<16 then MSB is the coarsest level 0
    return id >> (48 - (l+1)*3);
  }

  void assignIdx(size_t _i) { i = _i; }

  gridAtt3d() {}

  size_t id; // multi-level grid id

  size_t i; // index in the level

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

  sequence<gridVertex> P;

  sequence<size_t>** pointers;

  size_t L;

  floatT bSize;

  pt pMin;

  floatT maxSpan;

public:
  static constexpr size_t maxRange = 65535; // must match gridVertex.maxRange
  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  size_t numLevels() { return L; }

  gridVertex at(size_t l, size_t i) {
    if (l >= 0 && l < L) {
      return at(l+1, pointers[l]->at(i));
    } else if (l == L) {
      return P[i];
    } else {
      throw std::runtime_error("Invalid level for point access");
    }
  }

  sequence<size_t>* levelIdx(size_t l) {
    return pointers[l];
  }

  floatT boxSize(size_t l) {
    if (l >= 0 && l < L) {
      return bSize * pow(2, maxBit - l - 1);
    } else if (l == L) {
      return 0;
    } else {
      throw std::runtime_error("Invalid level for box size");
    }
  }

  floatT span() {
    return maxSpan;
  }

  pt getMin() {
    return pMin;
  }

  floatT levelSize(size_t l) {
    if (l >= 0 && l < L) {
      return pointers[l]->size() - 1;
    } else if (l == L) {
      return P.size();
    } else {
      throw std::runtime_error("Invalid level for level size");
    }
  }

  sequence<gridVertex> children(size_t l, size_t i) {
    if (l >= 0 && l < L) {
      //sequence<gridVertex> out(levelSize(l));
      size_t start = pointers[l]->at(i);
      size_t numChildren = pointers[l]->at(i+1) - start;
      sequence<gridVertex> out(numChildren);
      for (int i = 0; i < numChildren; ++i) {
	out[i] = at(l + 1, start + i);
	out[i].attribute.assignIdx(start + i); // Idx in the level
      }
      return out;
    } else if (l == L) {
      throw std::runtime_error("Invalid level L for children access");
    } else {
      throw std::runtime_error("Invalid level >L for children access");
    }
  }

  sequence<gridVertex> level(size_t l) {
    if (l >= 0 && l <= L) {
      sequence<gridVertex> out(levelSize(l));
      parallel_for(0, levelSize(l),
		   [&](size_t i){
				out[i] = at(l,i);
				out[i].attribute.assignIdx(i); // Idx in the level
			      });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  sequence<gridVertex> nextLevel(size_t l, sequence<size_t>& keep) {
    if (l >= 0 && l <= L && levelSize(l)+1 == keep.size()) {
      size_t ls = keep.size()-1;
      parallel_for(0, ls, [&](size_t i){
			    if (keep[i] > 0) {
			      keep[i] = pointers[l]->at(i+1) - pointers[l]->at(i);
			    } else
			      keep[i] = 0;
			  });
      size_t m = parlay::scan_inplace(keep.cut(0, ls));
      keep[ls] = m;

      sequence<gridVertex> out(m);
      parallel_for(0, ls, [&](size_t i){
			    if (keep[i] != keep[i+1]) {
			      size_t numPts = keep[i+1] - keep[i];
			      for (size_t j=0; j<numPts; ++j) {
				size_t nextLvlIdx = pointers[l]->at(i) + j;
				out[keep[i] + j] = at(l+1, nextLvlIdx);
				out[keep[i] + j].attribute.assignIdx(nextLvlIdx);
			      }
			    }
			  });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  ~gridTree3d() {
    for (size_t l=0; l<L; ++l)
      delete pointers[l];
    free(pointers);
  }

  gridTree3d(slice<pt*, pt*> _P) {
    L = maxBit;

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
    maxSpan = 1.01 * max(extrema[4]-extrema[5],
		  max((extrema[2]-extrema[3]),(extrema[0]-extrema[1])));
    size_t _maxRange = maxRange;
    bSize = maxSpan / maxRange;

    cout << "max-span = " << maxSpan << endl;
    // untranslated
    pMin[0] = extrema[1];
    pMin[1] = extrema[3];
    pMin[2] = extrema[5];

    cout << "grid-size = " << bSize << endl;

    P = sequence<gridVertex>(_P.size());
    parallel_for(0, P.size(), [&](size_t i){
				P[i] = gridVertex(_P[i].coords());
				P[i].attribute = gridAtt3d();
				P[i].attribute.computeId(_P[i], pMin, bSize);
				P[i].attribute.assignIdx(i);
			      });

    parlay::sort_inplace(make_slice(P), [&](gridVertex const& a, gridVertex const& b){
					  return a.attribute.id < b.attribute.id;});

    sequence<size_t> flag(P.size()+1);

    // Build the coarser L levels, from the finest to the coarsest (L-1 to 0)
    pointers = (sequence<size_t>**) malloc(sizeof(sequence<size_t>*)*L);
    for (int l=L-1; l>=0; --l) {
      flag[0] = 1;
      parallel_for(1, levelSize(l+1), [&](size_t i){
					if (at(l+1, i).attribute.gridId(l) !=
					    at(l+1, i-1).attribute.gridId(l)) {
					  flag[i] = 1;
					} else flag[i] = 0;
				      });
      size_t numGrids = parlay::scan_inplace(flag.cut(0, levelSize(l+1)));
      flag[levelSize(l+1)] = numGrids;
      pointers[l] = new sequence<size_t>(numGrids+1);
      parallel_for(0, levelSize(l+1), [&](size_t i){
				  if (flag[i] != flag[i+1]) {
				    pointers[l]->at(flag[i]) = i;
				  }
				});
      pointers[l]->at(numGrids) = levelSize(l+1);
      cout << "constructed level " << l << " of size " << numGrids << endl;
    }

  }
};

class gridOrigin;

template <class vertexT>
struct linkedFacet3d {
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

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
  seqT *keepList;

  void reassign(seqT *_seeList, gridOrigin* o) {
    delete seeList;
    delete keepList;

    //todo parallelize
    seeList = new seqT();
    keepList = new seqT();
    for (size_t i=0; i<_seeList->size(); ++i) {
      auto v = _seeList->at(i);
      if (o->visible(this, v) && a != v && b != v && c != v)
	seeList->push_back(v);
      else
	keepList->push_back(v);
    }

    delete _seeList;
  }

  void push_back(vertexT v, gridOrigin* o) {
    if (o->visible(this, v) && a != v && b != v && c != v)
      seeList->push_back(v);
    else
      keepList->push_back(v);
  }

  // Accesses the visible points (boxes)
  size_t numVisiblePts() { return seeList->size(); }
  inline vertexT& visiblePts(size_t i) { return seeList->at(i); }

  // Access both the visible and intersecting (boxes)
  inline size_t numPts() { return seeList->size() + keepList->size(); }
  vertexT& pts(size_t i) {
    if (i < seeList->size()) {
      return seeList->at(i);
    } else {
      return keepList->at(i - seeList->size());
    }
  }

  void clear() {
    seeList->clear();
    keepList->clear();
  }

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

    apex = parlay::max_element(seeList->cut(0, numVisiblePts()),
			       [&](vertexT aa, vertexT bb) {
				 return pargeo::signedVolumeX6(a, b, c, aa) <
				   pargeo::signedVolumeX6(a, b, c, bb);
			       });
    return apex;
  }

  linkedFacet3d(vertexT _a, vertexT _b, vertexT _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > _a.attribute.numericKnob)
      swap(b, c);

    seeList = new seqT();
    keepList = new seqT();

#ifdef RESERVE
    reservation = -1; // (unsigned)
#endif

    area = crossProduct3d(b-a, c-a);
  }

  ~linkedFacet3d() {
    delete seeList;
    delete keepList;
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
  using pt = pargeo::fpoint<3>;
  using floatT = typename pt::floatT;
  static constexpr floatT numericKnob = 1e-5;

  vertexT o; // translation origin

  vertexT pMin; // translated pmin of the grid

  floatT bSize; // gridSize of the level being processed

public:
  void setBoxSize(floatT _boxSize) {
    bSize = _boxSize;
  }

  template<class pt>
  void setMin(pt _pMin) {
    pMin = vertexT(_pMin.coords());
  }

  template<class pt>
  void setOrigin(pt _o) {
    o = vertexT(_o.coords());
    pMin = pMin - o; // translate the pmin by origin (same as other points in the hull)
  }

  inline vertexT get() {return o;};

  gridOrigin() {}

  template<class pt>
  gridOrigin(pt _o) {
    setOrigin(_o);
  }

  // Point visibility test
  inline bool visible(facetT* f, vertexT p) {
    return pargeo::signedVolumeX6(f->a, p, f->area) > numericKnob;
  }

  vertexT getCorner(vertexT p, vertexT myMin) {
    // todo compute from point grid id
    vertexT vt0 = p - pMin;
    vt0[0] = floor( vt0[0] / bSize ) * bSize;
    vt0[1] = floor( vt0[1] / bSize ) * bSize;
    vt0[2] = floor( vt0[2] / bSize ) * bSize;
    return vt0 + myMin;
  }

  // Input p needs to be untranslated
  inline void writeCorners(vertexT p, ofstream& os) {
    vertexT vt0 = getCorner(p, pMin + o);
    for (size_t i=0; i<8; ++i) {
      auto vt = vt0;
      vt[0] += (i & 1) * bSize;
      vt[1] += ((i>>1) & 1) * bSize;
      vt[2] += ((i>>2) & 1) * bSize;
      os << vt[0] << " " << vt[1] << " " << vt[2] << " " << i << endl;
    }
  }

  inline bool anyCornerVisible(facetT* f, vertexT p) {
    vertexT vt0 = getCorner(p, pMin);
    for (size_t i=0; i<8; ++i) {
      vertexT vt = vt0;
      vt[0] += (i & 1) * bSize;
      vt[1] += ((i>>1) & 1) * bSize;
      vt[2] += ((i>>2) & 1) * bSize;

      if (visible(f, vt))
	return true;
    }
    return false;
  }

  // This function returns true of any of the box corner of p is visible
  inline bool keep(facetT* f, vertexT p) {
    return anyCornerVisible(f, p);
  }
};
