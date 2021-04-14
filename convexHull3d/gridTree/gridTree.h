#pragma

#include <bitset>
#include "parlay/sequence.h"

using namespace parlay;

class gridAtt3d;

using gridVertex = pargeo::_point<3, pargeo::fpoint<3>::floatT, pargeo::fpoint<3>::floatT, gridAtt3d>;

struct gridAtt3d {
  using pt = pargeo::fpoint<3>;
  using floatT = pt::floatT;
  static constexpr typename pargeo::fpoint<3>::floatT numericKnob = 1e-5;

  gridAtt3d() {}

  size_t id; // grid id

  size_t i; // index, todo check if needed

  linkedFacet3d<gridVertex> *seeFacet;

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

template<class pt>
class gridTree3d {
  using floatT = gridVertex::floatT;
  static constexpr int dim = 3;
  static constexpr size_t maxRange = 65535;

  sequence<gridVertex> P;

  sequence<size_t>* pointers;

  size_t L;

  inline size_t shift(size_t x) {
    if (x > maxRange)
      throw std::runtime_error("Grid id exceed 16 bits, make min-grid size larger.");

    x = (x | (x << 16)) & 0x001f0000ff0000ff;
    x = (x | (x <<  8)) & 0x100f00f00f00f00f;
    x = (x | (x <<  4)) & 0x10c30c30c30c30c3;
    x = (x | (x <<  2)) & 0x1249249249249249;
    return x;
  }

  size_t computeId(pargeo::fpoint<3> p, pargeo::fpoint<3> pMin, floatT gSize) {
    size_t x = floor( (p[0] - pMin[0]) / gSize );
    size_t y = floor( (p[1] - pMin[1]) / gSize );
    size_t z = floor( (p[2] - pMin[2]) / gSize );

    x = shift(x);
    y = shift(y);
    z = shift(z);
    size_t id = x | (y << 1) | (z << 2);
    return id;
  }

  // Level 0 is the coarsest level
  size_t getLevel(gridVertex p, size_t l) {
    // the highest 16 bits are empty
    // starting taking 3-bit numbers
    return p.attribute.id >> (48 - (l+1)*3);
  }

public:

  sequence<size_t> levelIdx(size_t l) {
    return pointers[l];
  }

  sequence<gridVertex> level(size_t l) {
    sequence<gridVertex> out(pointers[l].size()-1);
    parallel_for(0, pointers[l].size()-1,
		 [&](size_t i){
		   out[i] = P[pointers[l][i]];
		 });
    return out;
  }

  gridVertex at(size_t l, size_t i) {
    return P[pointers[l][i]];
  }

  gridVertex at(size_t i) {
    return P[i];
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
    floatT gSize = max(extrema[4]-extrema[5],
		       max((extrema[2]-extrema[3]),(extrema[1]-extrema[0]))) / maxRange;

    pt pMin;
    pMin[0] = extrema[1];
    pMin[1] = extrema[3];
    pMin[2] = extrema[5];

    cout << "grid-size = " << gSize << endl;

    P = sequence<gridVertex>(_P.size());
    parallel_for(0, P.size(), [&](size_t i){
				P[i] = gridVertex(_P[i].coords());
				P[i].attribute = gridAtt3d();
				P[i].attribute.id = computeId(_P[i], pMin, gSize);
			      });

    parlay::sort_inplace(make_slice(P), [&](gridVertex const& a, gridVertex const& b){
					  return a.attribute.id < b.attribute.id;});

    sequence<size_t> flag(P.size()+1);

    // The the coarser L levels
    pointers = (sequence<size_t>*) malloc(sizeof(sequence<size_t>)*L);
    for (int l=0; l<L; ++l) {
      flag[0] = 1;
      parallel_for(1, P.size(), [&](size_t i){
				  if (getLevel(P[i], l) != getLevel(P[i-1], l)) {
				    flag[i] = 1;
				  } else flag[i] = 0;
				});
      flag[P.size()] = parlay::scan_inplace(flag.cut(0, flag.size()-1));

      pointers[l] = sequence<size_t>(flag[P.size()]+1);
      parallel_for(0, P.size(), [&](size_t i){
				  if (flag[i] != flag[i+1]) {
				    pointers[l][flag[i]] = i;
				  }
				});
      pointers[l][flag[P.size()]] = P.size();
      cout << "lvl " << l << " = " << flag[P.size()] << endl;
    }

  }
};
