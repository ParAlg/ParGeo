// This code is part of the project "ParGeo: A Library for Parallel Computational Geometry"
// Copyright (c) 2021-2022 Yiqiu Wang, Shangdi Yu, Laxman Dhulipala, Yan Gu, Julian Shun
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

#include "parlay/sequence.h"
#include "point.h"
#include <tuple>

namespace pargeo {

  template<int _dim, class _objT>
  class kdNode {

    using intT = int;
    using floatT = double;
    using pointT = _objT;
    using nodeT = kdNode<_dim, _objT>;

    intT id;

    int k;

    pointT pMin, pMax;

    nodeT* left;

    nodeT* right;

    nodeT* sib;

    parlay::slice<_objT**, _objT**> items;

    inline void minCoords(pointT& _pMin, pointT& p) {
      for(int i=0; i<_pMin.dim; ++i)
	_pMin[i] = std::min(_pMin[i], p[i]);
    }

    inline void maxCoords(pointT& _pMax, pointT& p) {
      for(int i=0; i<_pMax.dim; ++i)
	_pMax[i] = std::max(_pMax[i], p[i]);
    }

    void boundingBoxSerial();

    void boundingBoxParallel();

    intT splitItemSerial(floatT xM);

    inline bool itemInBox(pointT pMin1, pointT pMax1, _objT* item) {
      for(int i=0; i<_dim; ++i) {
	if (pMax1[i]<item->at(i) || pMin1[i]>item->at(i))
	  return false;
      }
      return true;
    }

    inline intT findWidest() {
      floatT xM = -1;
      for (int kk=0; kk<_dim; ++kk) {
	if (pMax[kk]-pMin[kk]>xM) {
	  xM = pMax[kk]-pMin[kk];
	  k = kk;}}
      return k;
    }

    void constructSerial(nodeT *space, intT leafSize);

    void constructParallel(nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize);

  public:
    using objT = _objT;

    static constexpr int dim  = _dim;

    inline nodeT* L() {return left;}

    inline nodeT* R() {return right;}

    inline nodeT* siblin() {return sib;}

    inline intT size() {return items.size();}

    inline _objT* operator[](intT i) {return items[i];}

    inline _objT* at(intT i) {return items[i];}

    inline bool isLeaf() {return !left;}

    inline _objT* getItem(intT i) {return items[i];}

    inline pointT getMax() {return pMax;}

    inline pointT getMin() {return pMin;}

    inline floatT getMax(int i) {return pMax[i];}

    inline floatT getMin(int i) {return pMin[i];}

    inline bool hasId() {
      return id != -1;
    }

    inline void setId(intT _id) {
      id = _id;
    }

    inline void resetId() {
      id = -1;
    }

    inline intT getId() {
      return id;
    }

    static const int boxInclude = 0;

    static const int boxOverlap = 1;

    static const int boxExclude = 2;

    inline floatT diag() {
      floatT result = 0;
      for (int d = 0; d < _dim; ++ d) {
	floatT tmp = pMax[d] - pMin[d];
	result += tmp * tmp;
      }
      return sqrt(result);
    }

    inline floatT lMax() {
      floatT myMax = 0;
      for (int d=0; d < _dim; ++d) {
	floatT thisMax = pMax[d] - pMin[d];
	if (thisMax > myMax) {
	  myMax = thisMax;}
      }
      return myMax;
    }

    inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2) {
      bool exclude = false;
      bool include = true; //1 include 2
      for(int i=0; i<_dim; ++i) {
	if (pMax1[i]<pMin2[i] || pMin1[i]>pMax2[i]) exclude = true;
	if (pMax1[i]<pMax2[i] || pMin1[i]>pMin2[i]) include = false;
      }
      if (exclude) return boxExclude;
      else if (include) return boxInclude;
      else return boxOverlap;
    }

    kdNode(parlay::slice<_objT**, _objT**> itemss, intT nn, nodeT *space,
	   parlay::slice<bool*, bool*> flags, intT leafSize=16);

    kdNode(parlay::slice<_objT**, _objT**> itemss, intT nn, nodeT *space,
	   intT leafSize=16);

  };

  // TODO move these two functions

  template <typename nodeT>
  double nodeDistance(nodeT* n1, nodeT* n2);

  template <typename nodeT>
  double nodeFarDistance(nodeT* n1, nodeT* n2);

  /* Kd-tree construction */

  template<int dim, class objT>
  kdNode<dim, objT>* buildKdTree(parlay::slice<objT*, objT*> P,
				 bool parallel = true,
				 size_t leafSize = 16,
				 parlay::sequence<objT*>* items = nullptr);

  template<int dim, class objT>
  kdNode<dim, objT>* buildKdTree(parlay::sequence<objT>& P,
				 bool parallel = true,
				 size_t leafSize = 16,
				 parlay::sequence<objT*>* items = nullptr);

  /* Kd-tree knn search */

  template<int dim, class objT>
  parlay::sequence<size_t> kdTreeKnn(parlay::sequence<objT> &queries,
				     size_t k,
				     kdNode<dim, objT>* tree = nullptr,
				     bool sorted = false);

  /* Kd-tree range search */

  template<int dim, typename objT>
  parlay::sequence<size_t> kdTreeRange(parlay::sequence<objT>& A,
				       kdNode<dim, objT>* tree,
				       objT query,
				       double radius);

  template<int dim, typename objT>
  parlay::sequence<size_t> kdTreeOrthRange(parlay::sequence<objT>& A,
					   kdNode<dim, objT>* tree,
					   objT query,
					   double halfLen);

  /* Bichromatic closest pair */

  template <typename nodeT>
  std::tuple<typename nodeT::objT*,
	     typename nodeT::objT*,
	     typename nodeT::objT::floatT> bccp(nodeT* n1, nodeT* n2);

  /* Well-separated pair decomposition */

  template <typename nodeT>
  struct wsp {
    nodeT* u;
    nodeT* v;
    wsp(nodeT* uu, nodeT* vv): u(uu), v(vv) {}
  };

  template<int dim>
  parlay::sequence<wsp<kdNode<dim, point<dim>>>>
  wspdParallel(kdNode<dim, point<dim>>* tree, double s = 2);

} // End namespace pargeo

#include "kdTreeImpl.h"
#include "kdTreeKnn.h"
#include "kdTreeRange.h"
#include "bccp.h"
#include "wspd.h"
