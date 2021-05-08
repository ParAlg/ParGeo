#pragma

#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace parlay;
using namespace pargeo;

//todo namespace
//todo figure out the perf gap vs no grid

template<typename pt>
class octTree {
  using floatT = typename pt::floatT;

  sequence<pt> P;

  sequence<size_t>** pointers;

  size_t L;

  floatT bSize; // the smallest box size

  sequence<floatT> boxSizes;

  sequence<size_t> levelIdx;

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

  size_t numLevels() { return L; }

  // read the representative point at level l offset i
  pt at(size_t l, size_t i) {
    if (l >= 0 && l < L) {
      return at(l+1, pointers[l]->at(i));
    } else if (l == L) {
      return P[i];
    } else {
      throw std::runtime_error("Invalid level for point access");
    }
  }

  // read the pointer value of level k, pointed by level l offset i
  size_t pointer(size_t l, size_t i, size_t k) {
    if (l >= 0 && l < k) {
      return pointer(l+1, pointers[l]->at(i), k);
    } else if (l == k) {
      return pointers[l]->at(i);
    } else {
      throw std::runtime_error("Invalid level for pointer access");
    }
  }

  size_t levelIndex(size_t l) {
    if (l >= 0 && l < L) {
      return levelIdx[l];
    } else if (l == L) {
      return maxBit;
    } else {
      throw std::runtime_error("Invalid level for box size");
    }
  }

  floatT boxSize(size_t l) {
    if (l >= 0 && l < L) {
      return boxSizes[l];
      //return bSize * pow(2, maxBit - l - 1);
    } else if (l == L) {
      return 0;
    } else {
      throw std::runtime_error("Invalid level for box size");
    }
    // if (l >= 0 && l < L) {
    //   return bSize * pow(2, maxBit - l - 1);
    // } else if (l == L) {
    //   return 0;
    // } else {
    //   throw std::runtime_error("Invalid level for box size");
    // }
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

  // return level l node i children
  // sequence<pt> children(size_t l, size_t i) {
  //   if (l >= 0 && l < L) {
  //     size_t start = pointers[l]->at(i);
  //     size_t numChildren = pointers[l]->at(i+1) - start;
  //     sequence<pt> out(numChildren);
  //     for (int i = 0; i < numChildren; ++i) {
  // 	out[i] = at(l + 1, start + i);
  // 	out[i].attribute.i = start + i; // Idx in the level
  //     }
  //     return out;
  //   } else if (l == L) {
  //     throw std::runtime_error("Invalid level L for children access");
  //   } else {
  //     throw std::runtime_error("Invalid level >L for children access");
  //   }
  // }

  sequence<pt> level(size_t l) {
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

  // Given level l mask, get level k points
  sequence<pt> getLevel(size_t l, sequence<size_t>& keep, size_t k) {
    if (l < k && levelSize(l)+1 == keep.size() && k >= 0 && k <= L) {
      size_t ls = keep.size()-1;
      parallel_for(0, ls, [&](size_t i){
			    if (keep[i] > 0)
			      keep[i] = pointer(l, i+1, k-1) - pointer(l, i, k-1);
			    else
			      keep[i] = 0;
			  });
      size_t m = parlay::scan_inplace(keep.cut(0, ls));
      keep[ls] = m;

      sequence<pt> out(m);
      parallel_for(0, ls, [&](size_t i){
			    if (keep[i] != keep[i+1]) {
			      size_t numPts = keep[i+1] - keep[i];
			      size_t base = pointer(l, i, k-1);
			      for (size_t j=0; j<numPts; ++j) {
				out[keep[i] + j] = at(k, base + j);
				out[keep[i] + j].attribute.i = base + j;
			      }
			    }
			  });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  sequence<pt> nextLevel(size_t l, sequence<size_t>& keep) {
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

      sequence<pt> out(m);
      parallel_for(0, ls, [&](size_t i){
			    if (keep[i] != keep[i+1]) {
			      size_t numPts = keep[i+1] - keep[i];
			      for (size_t j=0; j<numPts; ++j) {
				size_t nextLvlIdx = pointers[l]->at(i) + j;
				out[keep[i] + j] = at(l+1, nextLvlIdx);
				out[keep[i] + j].attribute.i = nextLvlIdx;
			      }
			    }
			  });
      return out;
    } else {
      throw std::runtime_error("Invalid level for level access");
    }
  }

  ~octTree() {
    for (size_t l=0; l<L; ++l)
      delete pointers[l];
    free(pointers);
  }

  template <typename pIn>
  void zorderSort(slice<pIn*, pIn*> _P) {
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

    P = sequence<pt>(_P.size());
    parallel_for(0, P.size(), [&](size_t i){
				P[i] = pt(_P[i].coords());
				P[i].attribute = gridAtt3d();
				P[i].attribute.id = computeId(_P[i], pMin, bSize);
				P[i].attribute.i = i;
			      });

    parlay::sort_inplace(make_slice(P), [&](pt const& a, pt const& b){
					  return a.attribute.id < b.attribute.id;});
  }

  template <typename pIn>
  octTree(slice<pIn*, pIn*> _P) {
    zorderSort(_P);
    constructAllLevels();
  }

  size_t constructAllLevels() {
    sequence<size_t> flag(P.size()+1);
    // Build the coarser maxBit levels, from the finest to the coarsest (maxBit-1 to 0)

    L = maxBit;
    pointers = (sequence<size_t>**) malloc(sizeof(sequence<size_t>*)*maxBit);
    boxSizes = sequence<floatT>(maxBit);
    levelIdx = sequence<size_t>(maxBit);

    floatT s = bSize;

    for (int l = maxBit-1; l >= 0; -- l) {
      timer t; t.start();

      boxSizes[l] = s;
      s *= 2;
      levelIdx[l] = l;

      // access indirections here
      flag[0] = 1;
      parallel_for(1, levelSize(l+1), [&](size_t i){
					if (gridId(l, at(l+1, i).attribute.id) !=
					    gridId(l, at(l+1, i-1).attribute.id)) {
					  flag[i] = 1;
					} else flag[i] = 0;
				      });
      size_t numGrids = parlay::scan_inplace(flag.cut(0, levelSize(l+1)));
      flag[levelSize(l+1)] = numGrids;
      cout << " scanning-time = " << t.get_next() << endl;

      pointers[l] = new sequence<size_t>(numGrids+1);
      parallel_for(0, levelSize(l+1), [&](size_t i){
					if (flag[i] != flag[i+1]) {
					  pointers[l]->at(flag[i]) = i;
					}
				      });
      pointers[l]->at(numGrids) = levelSize(l+1);
      cout << " pointing-time =  " << t.stop() << endl;
      cout << "num-grids = " << numGrids << endl;
    }
    return maxBit;
  }

  // toConstruct: the level numbers in the reverse order
  template <typename pIn>
  octTree(slice<pIn*, pIn*> _P, sequence<size_t>& toConstruct) {
    zorderSort(_P);
    constructLevels(toConstruct);
  }

  size_t constructLevels(sequence<size_t>& toConstruct) {
    sequence<size_t> flag(P.size()+1);
    // Build the coarser maxBit levels, from the finest to the coarsest (maxBit-1 to 0)

    L = toConstruct.size();

    pointers = (sequence<size_t>**) malloc(sizeof(sequence<size_t>*)*L);
    boxSizes = sequence<floatT>(L);
    levelIdx = sequence<size_t>(L);

    floatT s = bSize;

    auto construct = [&](int l) {
		       for (auto x: toConstruct) {
			 if (x == l) return true;
		       }
		       return false;
		     };

    size_t lIdx = toConstruct.size() - 1;

    for (int l = maxBit-1; l >= 0; -- l) {
      cout << l << endl;

      if (!construct(l)) {
	s *= 2;
	continue;
      }

      timer t; t.start();

      boxSizes[lIdx] = s;
      cout << "level " << l << ": " << s << endl;
      s *= 2;
      levelIdx[lIdx] = l;

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
