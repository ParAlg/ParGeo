#pragma

#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace parlay;
using namespace pargeo;

//todo namespace

template<typename pt>
class octTree {
  using floatT = typename pt::floatT;

  sequence<size_t>** pointers;

  sequence<pt> P;

  size_t L;

  sequence<floatT> boxSizes;

  pt pMin;

  floatT maxSpan;

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

  template <typename pt1, typename pt2>
  size_t computeId(pt1 p, pt2 pMin, floatT boxSize) {
    size_t x = floor( (p[0] - pMin[0]) / boxSize );
    size_t y = floor( (p[1] - pMin[1]) / boxSize );
    size_t z = floor( (p[2] - pMin[2]) / boxSize );

    x = shift(x);
    y = shift(y);
    z = shift(z);
    return x | (y << 1) | (z << 2);
  }

  inline size_t gridId(size_t l, size_t id) {
    // The highest 16 bits are empty, each level is 3-bits
    // <<16 then MSB is the coarsest level 0
    return id >> (48 - (l+1)*3);
  }

public:

  static constexpr size_t maxRange = 65535;

  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  // Constructor
  // - toConstruct: the levels to construct in decreasing id
  template <typename pIn>
  octTree(slice<pIn*, pIn*> _P, sequence<size_t>& toConstruct) {
    floatT bSize = zorderSort(_P);
    constructLevels(toConstruct, bSize);
  }

  ~octTree() {
    for (size_t l=0; l<L; ++l)
      delete pointers[l];
    free(pointers);
  }

  size_t numLevels() { return L; }

  floatT boxSize(size_t l) {
    if (l >= 0 && l < L) {
      return boxSizes[l];
    } else if (l == L) {
      return 0;
    } else {
      throw std::runtime_error("Invalid level for box size");
    }
  }

  pt coordinateMin() {
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

  // Read representative point, level l offset i
  pt at(size_t l, size_t i) {
    if (l >= 0 && l < L) {
      return at(l+1, pointers[l]->at(i));
    } else if (l == L) {
      return P[i];
    } else {
      throw std::runtime_error("Invalid level for point access");
    }
  }

  sequence<pt> levelSample(size_t l) {
    if (l >= 0 && l <= L) {
      sequence<pt> out(levelSize(l));
      parallel_for(0, levelSize(l),
		   [&](size_t i){
				out[i] = at(l,i);
				out[i].attribute.i = i; // Idx in the level
			      });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  template <typename ptOut>
  sequence<ptOut> getPoints(size_t l, sequence<size_t>& mask) {
    if (levelSize(l)+1 == mask.size() && l <= L) {
      size_t ls = mask.size()-1;

      parallel_for(0, ls-1, [&](size_t i) {
			    if (mask[i] > 0) {
			      mask[i] = at(l,i+1).attribute.i - at(l,i).attribute.i;
			    } else {
			      mask[i] = 0;}
			  });
      mask[ls-1] = P.size() - at(l,ls-1).attribute.i;

      size_t m = parlay::scan_inplace(mask.cut(0, ls));
      mask[ls] = m;

      sequence<ptOut> out(m);
      parallel_for(0, ls, [&](size_t i){
			    if (mask[i] != mask[i+1]) {
			      size_t numPts = mask[i+1] - mask[i];
			      size_t base = at(l,i).attribute.i;
			      for (size_t j=0; j<numPts; ++j) {
				out[mask[i] + j] = ptOut(P[base + j].coords());
			      }
			    }
			  });
      return out;
    } else {
      throw std::runtime_error("Invalid level for base point access");
    }
  }

  sequence<pt> nextLevelSample(size_t l, sequence<size_t>& mask) {
    if (l >= 0 && l <= L && levelSize(l)+1 == mask.size()) {
      size_t ls = mask.size()-1;
      parallel_for(0, ls, [&](size_t i){
			    if (mask[i] > 0) {
			      mask[i] = pointers[l]->at(i+1) - pointers[l]->at(i);
			    } else
			      mask[i] = 0;
			  });
      size_t m = parlay::scan_inplace(mask.cut(0, ls));
      mask[ls] = m;

      sequence<pt> out(m);
      parallel_for(0, ls, [&](size_t i){
			    if (mask[i] != mask[i+1]) {
			      size_t numPts = mask[i+1] - mask[i];
			      for (size_t j=0; j<numPts; ++j) {
				size_t nextLvlIdx = pointers[l]->at(i) + j;
				out[mask[i] + j] = at(l+1, nextLvlIdx);
				out[mask[i] + j].attribute.i = nextLvlIdx;
			      }
			    }
			  });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  template <typename pIn>
  floatT zorderSort(slice<pIn*, pIn*> _P) {
    timer t; t.start();
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
    floatT bSize = maxSpan / maxRange;

    // untranslated
    pMin[0] = extrema[1];
    pMin[1] = extrema[3];
    pMin[2] = extrema[5];

    P = sequence<pt>(_P.size());
    parallel_for(0, P.size(), [&](size_t i){
				P[i] = pt(_P[i].coords());
				P[i].attribute = gridAtt3d();
				P[i].attribute.id = computeId(_P[i], pMin, bSize);
			      });
    cout << "grid-init-time = " << t.get_next() << endl;

    parlay::sort_inplace(make_slice(P), [&](pt const& a, pt const& b){
					  return a.attribute.id < b.attribute.id;});
    parallel_for(0, P.size(), [&](size_t i){
				P[i].attribute.i = i;
			      });
    cout << "grid-sort-time = " << t.stop() << endl;
    return bSize;
  }

  size_t constructLevels(sequence<size_t>& toConstruct, floatT bSize) {
    { // check input
      size_t l = maxBit;
      for (auto x: toConstruct) {
	if (x < 0 || x > maxBit-1 || x >= l)
	  throw std::runtime_error("level # must be [0,15], in decreasing order");
	l = x;
      }
    }

    sequence<size_t> flag(P.size()+1);
    // Build the coarser maxBit levels, from the finest to the coarsest (maxBit-1 to 0)

    L = toConstruct.size();

    pointers = (sequence<size_t>**) malloc(sizeof(sequence<size_t>*)*L);
    boxSizes = sequence<floatT>(L);

    floatT s = bSize;

    auto construct = [&](int l) {
		       for (auto x: toConstruct) {
			 if (x == l) return true;
		       }
		       return false;
		     };

    size_t lIdx = toConstruct.size() - 1;

    for (int l = maxBit-1; l >= 0; -- l) {

      if (!construct(l)) {
	s *= 2;
	continue;
      }

      timer t; t.start();

      boxSizes[lIdx] = s;
      cout << "--- level " << l << ": " << s << endl;
      s *= 2;

      flag[0] = 1;
      parallel_for(1, levelSize(lIdx+1), [&](size_t i){
					if (gridId(l, at(lIdx+1, i).attribute.id) !=
					    gridId(l, at(lIdx+1, i-1).attribute.id)) {
					  flag[i] = 1;
					} else flag[i] = 0;
				      });

      size_t numGrids = parlay::scan_inplace(flag.cut(0, levelSize(lIdx+1)));

      flag[levelSize(lIdx+1)] = numGrids;
      cout << " scanning-time = " << t.get_next() << endl;

      pointers[lIdx] = new sequence<size_t>(numGrids+1);
      parallel_for(0, levelSize(lIdx+1), [&](size_t i){
					if (flag[i] != flag[i+1]) {
					  pointers[lIdx]->at(flag[i]) = i;
					}
				      });
      pointers[lIdx]->at(numGrids) = levelSize(lIdx+1);
      cout << " pointing-time =  " << t.stop() << endl;
      cout << " num-grids = " << numGrids << endl;

      lIdx --;
    }
    return L;
  }

};
