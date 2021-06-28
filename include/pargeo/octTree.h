#pragma once

#include <math.h>
#include <tuple>
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include <bitset>

namespace pargeo {

  template<class objT, class pointT>
  class octTree {

  protected:
    using floatT = typename pointT::floatT;
    static constexpr size_t leafThresh = 16;
    static constexpr int dim = 3;

    class node {
    public:
      size_t s, e;

      node(): s(-1), e(-1) {}

      inline size_t size() {return e - s;}

      inline bool isEmpty() {return s < 0;}

      inline bool isLeaf() {
	return !isEmpty() && size() <= leafThresh;
      }
    };

    class pointer {
    public:
      objT p;
      size_t id;
      pointer(): id(-1) {}
      pointer(objT _p, size_t _id):
	p(_p), id(_id) {}
      bool isEmpty() {return id < 0;}
    };

    parlay::sequence<pointer> items;

    pointT pMin;

    floatT maxSpan, bSize;

    node** nodes;
    floatT* boxSizes;

    size_t L, maxNodes, wasteBits;

    inline bool lastLevel(size_t l) {
      return l == L;
    }

    inline size_t shift(size_t x) {
      x = (x | (x << 16)) & 0x001f0000ff0000ff;
      x = (x | (x <<  8)) & 0x100f00f00f00f00f;
      x = (x | (x <<  4)) & 0x10c30c30c30c30c3;
      x = (x | (x <<  2)) & 0x1249249249249249;
      return x;
    }

    template <class pt1>
    size_t id(pt1 p) {
      size_t x = floor( (p[0] - pMin[0]) / bSize );
      size_t y = floor( (p[1] - pMin[1]) / bSize );
      size_t z = floor( (p[2] - pMin[2]) / bSize );

      x = shift(x);
      y = shift(y);
      z = shift(z);
      return x | (y << 1) | (z << 2);
    }

    inline size_t octantId(size_t l, size_t id) {
      return (id >> ((L - 1 - l) * dim + wasteBits)) & size_t(7);
    }

    void initNodes() {
      maxNodes = 0;
      nodes = (node**) malloc(sizeof(node*) * (L + 1));
      for (size_t l = 0; l < L + 1; ++ l) { // todo OTF
        size_t s = pow(double(2), double(dim * l));
	nodes[l] = (node*) malloc(sizeof(node) * s);
	parlay::parallel_for(0, s, [&](size_t i) {
				     nodes[l][i] = node();
				   });
	maxNodes += s;
      }

      wasteBits = dim * (16 - L);

      if (L > 9) // upper limit is maxBits
	throw std::runtime_error("octtree is too deep, decrease L");
    }

    void init(parlay::slice<objT*, objT*> P) {
      initNodes();

      std::atomic<floatT> extrema[6];

      for (int i=0; i<3; ++i) {
	extrema[i*2] = P[0][i];
	extrema[i*2+1] = P[0][i];
      }

      parlay::parallel_for(0, P.size(), [&](size_t i){
          parlay::write_max(&extrema[0], P[i][0], std::less<floatT>());
	  parlay::write_min(&extrema[1], P[i][0], std::less<floatT>());
	  parlay::write_max(&extrema[2], P[i][1], std::less<floatT>());
	  parlay::write_min(&extrema[3], P[i][1], std::less<floatT>());
	  parlay::write_max(&extrema[4], P[i][2], std::less<floatT>());
	  parlay::write_min(&extrema[5], P[i][2], std::less<floatT>());
	});

      maxSpan = 1.01 * std::max(extrema[4]-extrema[5],
				std::max((extrema[2]-extrema[3]),
					 (extrema[0]-extrema[1])));

      floatT tmp = maxSpan;
      boxSizes = (floatT*) malloc(sizeof(floatT) * (L + 1));

      bSize = maxSpan / maxRange;

      for (size_t i = 0; i < L + 1; ++ i) {
	boxSizes[i] = tmp;
	tmp /= 2;
      }

      pMin[0] = extrema[1];
      pMin[1] = extrema[3];
      pMin[2] = extrema[5];
      // std::cout << "pMin = " << pMin << "\n";

      items = std::move(parlay::tabulate(P.size(),
					 [&](size_t i) {
					   return pointer(P[i], id(P[i]));
					 }));

    }

    inline size_t writeAdd(std::atomic<size_t>* a, size_t b) {
      size_t newV, oldV;
      do {
	oldV = a->load();
	newV = oldV + b;
      } while (!std::atomic_compare_exchange_weak(a, &oldV, newV));
      return oldV;
    }

    /* r is initialized before calling */
    void construct(size_t l, size_t nodeIdx, node *r) {
      // std::cout << "\nconstruct " << l << ": " << r->s << " -- " << r->e << "\n";

      if (r->isLeaf()) {
      	return;
      }

      auto Q = items.cut(r->s, r->e);

      parlay::integer_sort_inplace(Q, [&](pointer p) {
					return octantId(l, p.id);
				      });

      // size_t prev = octantId(l, Q[0].id);
      // bool error = false;
      // for (auto q: Q) {
      // 	if (prev > octantId(l, q.id)) error = true;
      // 	std::cout << octantId(l, q.id) << " ";
      // }
      // std::cout << "\n";
      // if (error) abort();

      size_t offsets[9]; // todo 8
      for (size_t i = 0; i < 9; ++ i) offsets[i] = 0;
      parlay::parallel_for(1, Q.size(), [&](size_t i) {
				  if (octantId(l, Q[i].id) !=
				      octantId(l, Q[i - 1].id)) {
				    offsets[octantId(l, Q[i].id)] = i;
				  }
				});
      for (size_t i = 1; i < 8; ++ i) {
	if (offsets[i] == 0) offsets[i] = offsets[i - 1];
      }

      // std::cout << "offsets = ";
      for (size_t i = 0; i < 8; ++ i) {
	offsets[i] += r->s;
	// std::cout << offsets[i] << ", ";
      }
      offsets[8] = r->e;
      // std::cout << "\n";

      if (lastLevel(l)) {
	// std::cout << "end construction\n";
	return;
      } else {

	if (r->e - r->s > 2000) {

	  parlay::parallel_for(0, 8, [&](size_t i){
			       if (offsets[i + 1] - offsets[i] > 0) {
				 size_t nodeIdxNext = nodeIdx << dim;
				 nodeIdxNext = nodeIdxNext + i;
				 node *q = &nodes[l + 1][nodeIdxNext];
				 q->s = offsets[i];
				 q->e = offsets[i + 1];
				 construct(l + 1, nodeIdxNext, q);
			       }			       
			     }, 1);

	} else {

	  for (size_t i = 0; i < 8; ++ i) {
	    if (offsets[i + 1] - offsets[i] > 0) {
	      size_t nodeIdxNext = nodeIdx << dim;
	      nodeIdxNext = nodeIdxNext + i;
	      node *q = &nodes[l + 1][nodeIdxNext];
	      q->s = offsets[i];
	      q->e = offsets[i + 1];
	      construct(l + 1, nodeIdxNext, q);
	    }
	  }

	}

      };
    }

  public:
    static constexpr size_t maxRange = 65535;
    static constexpr size_t maxBits = 16; // 2**maxBits = maxRange

    void printBits(size_t y) {
      std::bitset<64> x(y);
      std::cout << x << '\n';
    }

    std::tuple<objT, objT> getBox(size_t l, size_t nodeIdx) {
      // todo recover directly from nodeIdx
      node* r = &nodes[l][nodeIdx];
      objT p = items[r->s].p;
      floatT side = boxSize(l);
      size_t x = floor( (p[0] - pMin[0]) / side );
      size_t y = floor( (p[1] - pMin[1]) / side );
      size_t z = floor( (p[2] - pMin[2]) / side );

      objT bMin;
      bMin[0] = x * side;
      bMin[1] = y * side;
      bMin[2] = z * side;
      bMin = bMin + pMin;
      objT bMax;
      bMax[0] = bMin[0] + side;
      bMax[1] = bMin[1] + side;
      bMax[2] = bMin[2] + side;
      return std::tuple(bMin, bMax);
    }

    inline floatT boxSize(size_t l) {
      return boxSizes[l];
    }

    void info() {
      std::cout << "oct-tree:\n";
      std::cout << " max #-nodes = " << maxNodes << "\n";
      std::cout << " mem-size = " << (double(maxNodes) * sizeof(node) / 1000000) << " MB\n";
      std::cout << " bit-used = " << dim * L << "/64\n";
    }

    octTree(parlay::slice<objT*, objT*> P, size_t _L): L(_L) {
      init(P);
      // info();
      nodes[0][0].s = 0;
      nodes[0][0].e = P.size();
      construct(0, 0, &nodes[0][0]);
    }

    ~octTree() {
      for (size_t l = 0; l < L + 1; ++ l) free(nodes[l]);
      free(nodes);
    }
  };

} // End namespace pargeo
