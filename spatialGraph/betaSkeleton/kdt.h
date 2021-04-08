// This code is part of the "Pargeo" project
// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <limits> // numeric_limits
//#include <algorithm> // nth_element
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "geometry/point.h"

namespace skeletonKdt {

  using namespace std;

  template<int dim, class objT>
  class kdNode {
    typedef int intT;
    typedef double floatT;
    typedef pargeo::point<dim> pointT;
    typedef kdNode<dim, objT> nodeT;

    static const int boxInclude = 0;
    static const int boxOverlap = 1;
    static const int boxExclude = 2;

    // Data fields

    int k;
    pointT pMin, pMax;
    parlay::slice<objT**, objT**> items;
    intT n;
    nodeT* left;
    nodeT* right;
    nodeT* sib;

    // Methods

    inline void minCoords(pointT& _pMin, pointT& p) {
      for(int i=0; i<_pMin.dim; ++i)
	_pMin[i] = min(_pMin[i], p[i]);
    }

    inline void maxCoords(pointT& _pMax, pointT& p) {
      for(int i=0; i<_pMax.dim; ++i)
	_pMax[i] = max(_pMax[i], p[i]);
    }

    inline void boundingBoxSerial() {
      pMin = pointT(items[0]->coords());
      pMax = pointT(items[0]->coords());
      for(intT i=0; i<n; ++i) {
	minCoords(pMin, items[i][0]);
	maxCoords(pMax, items[i][0]);
      }
    }

    inline void boundingBoxParallel() {
      intT P = parlay::num_workers()*8;
      intT blockSize = (n+P-1)/P;
      pointT localMin[P];
      pointT localMax[P];
      for (intT i=0; i<P; ++i) {
	localMin[i] = pointT(items[0]->coords());
	localMax[i] = pointT(items[0]->coords());}
      parlay::parallel_for(0, P,
			   [&](intT p) {
			     intT s = p*blockSize;
			     intT e = min((intT)(p+1)*blockSize,n);
			     for (intT j=s; j<e; ++j) {
			       minCoords(localMin[p], items[j][0]);
			       maxCoords(localMax[p], items[j][0]);}
			   });
      pMin = pointT(items[0]->coords());
      pMax = pointT(items[0]->coords());
      for(intT p=0; p<P; ++p) {
	minCoords(pMin, localMin[p]);
	maxCoords(pMax, localMax[p]);}
    }

    inline intT splitItemSerial(floatT xM) {
      if (n < 2) {
	throw std::runtime_error("Error, kdTree splitting singleton.");}
      intT lPt = 0;
      intT rPt = n-1;
      while (lPt < rPt) {
	if (items[lPt]->at(k)>=xM) {
	  while (items[rPt]->at(k)>=xM && lPt < rPt) {
	    rPt--;
	  }
	  if (lPt < rPt) {
	    swap(items[lPt], items[rPt]);
	    rPt--; }
	  else { break;}
	}
	lPt++;
      }
      if (items[lPt]->at(k) < xM) lPt++;
      return lPt;
    }

    intT findWidest() {
      floatT xM = -1;
      for (int kk=0; kk<dim; ++kk) {
	if (pMax[kk]-pMin[kk]>xM) {
	  xM = pMax[kk]-pMin[kk];
	  k = kk;}}
      return k;
    }

    void constructSerial(nodeT *space, intT leafSize) {
      boundingBoxSerial();
      sib = NULL;
      if (n <= leafSize) {
	left = NULL; right = NULL;
      } else {
	intT k = findWidest();
	floatT xM = (pMax[k]+pMin[k])/2;

	// Split items by xM (serial)
	intT median = splitItemSerial(xM);

	if (median == 0 || median == n) {median = ceil(n/2.0);}

	if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) {
	  throw std::runtime_error("Error, kdNode overwrite.");
	}

	// Recursive construction
	space[0] = nodeT(items.cut(0, median), median, space+1, leafSize);
	space[2*median-1] = nodeT(items.cut(median, n), n-median, space+2*median, leafSize);
	left = space;
	right = space+2*median-1;
	left->sib = right;
	right->sib = left;
      }
    }

    void constructParallel(nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize) {
      boundingBoxParallel();

      sib = NULL;
      if (n <= leafSize) {
	left = NULL; right = NULL;
      } else {
	intT k = findWidest();
	floatT xM = (pMax[k]+pMin[k])/2;

	// Split items by xM in dim k (parallel)
	parlay::parallel_for(0, n,
			     [&](intT i) {
			       if (items[i]->at(k)<xM) flags[i]=1;
			       else flags[i] = 0;});
	auto mySplit = parlay::internal::split_two(items, flags);
	auto splited = mySplit.first;
	intT median = mySplit.second;
	parlay::parallel_for(0, n, [&](intT i) {items[i] = splited[i];}); // Copy back

	if (median == 0 || median == n) {median = (n/2.0);}

	if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) {
	  throw std::runtime_error("Error, kdNode overwrite.");
	}

	// Recursive construction
	parlay::par_do([&](){space[0] = nodeT(items.cut(0, median), median, space+1, flags.cut(0, median), leafSize);},
		       [&](){space[2*median-1] = nodeT(items.cut(median, n), n-median, space+2*median, flags.cut(median, n), leafSize);});
	left = space;
	right = space+2*median-1;
	left->sib = right;
	right->sib = left;
      }
    }

  public:

    inline nodeT* L() {return left;}

    inline nodeT* R() {return right;}

    inline nodeT* siblin() {return sib;}//todo

    inline intT size() {return n;}

    inline objT* operator[](intT i) {return items[i];}

    inline void setEmpty() {n=-1;}

    inline bool isEmpty() {return n<0;}

    inline bool isLeaf() {return !left;}//check

    inline objT* getItem(intT i) {return items[i];}

    inline pointT getMax() {return pMax;}

    inline pointT getMin() {return pMin;}

    kdNode(parlay::slice<objT**, objT**> itemss, intT nn, nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize=16): items(itemss), n(nn) {
      if (n>2000) constructParallel(space, flags, leafSize);
      else constructSerial(space, leafSize);
    }

    kdNode(parlay::slice<objT**, objT**> itemss, intT nn, nodeT *space, intT leafSize=16): items(itemss), n(nn) {
      constructSerial(space, leafSize);//todo get rid of intT n
    }

    bool nonEmptyLuneHelper(pointT cMin1, pointT cMax1, pointT cMin2, pointT cMax2, objT& p1, double r1, objT& p2, double r2, objT* e1, objT* e2);
    bool nonEmptyLune(objT& p1, double r1, objT& p2, double r2, objT* e1, objT* e2);
  };

  template<int dim, class objT>
  kdNode<dim, objT>* buildKdt(parlay::sequence<objT>& P, bool parallel=true, bool noCoarsen=false) {
    typedef kdNode<dim, objT> nodeT;

    size_t n = P.size();

    auto items = parlay::sequence<objT*>(n);
    parlay::parallel_for(0, n, [&](size_t i) {items[i]=&P[i];});
    parlay::slice<objT**, objT**> itemSlice = parlay::slice(items.begin(), items.end());

    auto root = (nodeT*) malloc(sizeof(nodeT)*(2*n-1));
    parlay::parallel_for(0, 2*n-1, [&](size_t i) {
				     root[i].setEmpty();
				   });

    if (parallel) {
      auto flags = parlay::sequence<bool>(n);
      auto flagSlice = parlay::slice(flags.begin(), flags.end());
      root[0] = nodeT(itemSlice, n, root+1, flagSlice, noCoarsen ? 1 : 16);
    } else {
      root[0] = nodeT(itemSlice, n, root+1, noCoarsen ? 1 : 16);
    }

    return root;
  }

  /////////////////////////////////////////////////////////
  // Routines specialized for beta skeleton searches in 2D
  /////////////////////////////////////////////////////////

  /*
    - Circle 1 bounding box: cMin1, cMax1
    - Circle 2 bounding box: cMin2, cMax2
    - Circle 1: p1, r1
    - Circle 2: p2, r2
    - Edge vertices: e1, e2
   */
  template<int dim, class objT>
  bool kdNode<dim, objT>::nonEmptyLuneHelper(pointT cMin1, pointT cMax1, pointT cMin2, pointT cMax2, objT& p1, double r1, objT& p2, double r2, objT* e1, objT* e2) {
    enum stateT {Include, Disjoint, Intersect};

    auto boxCompare = [&](pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2) {
      bool exclude = false;
      bool include = true;//1 include 2
      for(int i=0; i<dim; ++i) {
	if (pMax1[i]<pMin2[i] || pMin1[i]>pMax2[i]) exclude = true;
	if (pMax1[i]<pMax2[i] || pMin1[i]>pMin2[i]) include = false;
      }
      if (exclude) return Disjoint;
      else if (include) return Include;
      else return Intersect;
    };

    stateT relation1 = boxCompare(cMin1, cMax1, getMin(), getMax());
    stateT relation2 = boxCompare(cMin2, cMax2, getMin(), getMax());

    if(relation1 == Disjoint || relation2 == Disjoint) {
      return false;
    /* } else if (relation1 == Include && relation2 == Include) { */
    /*   return true; */
    } else {
      // The box at least intersects with both circles
      if (isLeaf()) {
	for (intT i=0; i < size(); ++ i) {
	  objT* p = getItem(i);
	  if (p == e1 || p == e2) continue; // Disregard edge end points
	  float dist1 = p1.dist(*p);
	  float dist2 = p2.dist(*p);
	  if (dist1 <= r1 && dist2 <= r2)
	    return true;
	}
	return false;
      } else {
	return
	  L()->kdNode<dim, objT>::nonEmptyLuneHelper(cMin1, cMax1, cMin2, cMax2, p1, r1, p2, r2, e1, e2) ||
	  R()->kdNode<dim, objT>::nonEmptyLuneHelper(cMin1, cMax1, cMin2, cMax2, p1, r1, p2, r2, e1, e2);
      }
    }
  }

  /* Returns whether the lune defined by circle 1 and circle is non-empty,
     not considering he edge vertices e1 and e2
    - Circle 1: p1, r1
    - Circle 2: p2, r2
    - Edge vertices: e1, e2
  */
  template<int dim, class objT>
  bool kdNode<dim, objT>::nonEmptyLune(objT& p1, double r1, objT& p2, double r2, objT* e1, objT* e2) {
    pointT cMin1, cMax1, cMin2, cMax2;
    cMin1[0] = p1[0] - r1; cMin1[1] = p1[1] - r1;
    cMax1[0] = p1[0] + r1; cMax1[1] = p1[1] + r1;
    cMin2[0] = p2[0] - r2; cMin2[1] = p2[1] - r2;
    cMax2[0] = p2[0] + r2; cMax2[1] = p2[1] + r2;
    return nonEmptyLuneHelper(cMin1, cMax1, cMin2, cMax2, p1, r1, p2, r2, e1, e2);
  }

} // End namespace
