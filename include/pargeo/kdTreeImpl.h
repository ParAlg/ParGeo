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

#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "kdTree.h"

namespace pargeo {

  namespace kdTreeInternal {

    template <typename In_Seq, typename Bool_Seq>
    auto split_two(In_Seq const &In,
		   Bool_Seq const &Fl,
		   parlay::flags fl = parlay::no_flag)
      -> std::pair<parlay::sequence<typename In_Seq::value_type>, size_t> {

      using namespace parlay;
      using namespace parlay::internal;

      using T = typename In_Seq::value_type;
      size_t n = In.size();
      size_t l = num_blocks(n, _block_size);
      sequence<size_t> Sums(l);
      sliced_for(n, _block_size,
		 [&](size_t i, size_t s, size_t e) {
		   size_t c = 0;
		   for (size_t j = s; j < e; j++) c += (Fl[j] == false);
		   Sums[i] = c;
		 },
		 fl);
      size_t m = scan_inplace(make_slice(Sums), addm<size_t>());
      sequence<T> Out = sequence<T>::uninitialized(n);
      sliced_for(n, _block_size,
		 [&](size_t i, size_t s, size_t e) {
		   size_t c0 = Sums[i];
		   size_t c1 = s + (m - c0);
		   for (size_t j = s; j < e; j++) {
		     if (Fl[j] == false)
		       assign_uninitialized(Out[c0++], In[j]);
		     else
		       assign_uninitialized(Out[c1++], In[j]);
		   }
		 },
		 fl);
      return std::make_pair(std::move(Out), m);
    }
  } // End namespace kdTreeInternal

  template<int _dim, class _objT>
  void kdNode<_dim, _objT>::boundingBoxSerial() {
    pMin = pointT(items[0]->coords());
    pMax = pointT(items[0]->coords());
    for(intT i=0; i<size(); ++i) {
      minCoords(pMin, items[i][0]);
      maxCoords(pMax, items[i][0]);
    }
  }

  template<int _dim, class _objT>
  void kdNode<_dim, _objT>::boundingBoxParallel() {
    intT P = parlay::num_workers()*8;
    intT blockSize = (size()+P-1)/P;
    pointT localMin[P];
    pointT localMax[P];
    for (intT i=0; i<P; ++i) {
      localMin[i] = pointT(items[0]->coords());
      localMax[i] = pointT(items[0]->coords());}
    parlay::parallel_for(0, P,
			 [&](intT p) {
			   intT s = p*blockSize;
			   intT e = std::min((intT)(p+1)*blockSize, size());
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

  template<int _dim, class _objT>
  typename kdNode<_dim, _objT>::intT
  kdNode<_dim, _objT>::splitItemSerial(floatT xM) {
    if (size() < 2) {
      throw std::runtime_error("Error, kdTree splitting singleton.");}
    intT lPt = 0;
    intT rPt = size()-1;
    while (lPt < rPt) {
      if (items[lPt]->at(k)>=xM) {
	while (items[rPt]->at(k)>=xM && lPt < rPt) {
	  rPt--;
	}
	if (lPt < rPt) {
	  std::swap(items[lPt], items[rPt]);
	  rPt--; }
	else { break;}
      }
      lPt++;
    }
    if (items[lPt]->at(k) < xM) lPt++;
    return lPt;
  }

  template<int _dim, class _objT>
  void kdNode<_dim, _objT>::constructSerial(nodeT *space, intT leafSize) {
    boundingBoxSerial();
    sib = NULL;
    if (size() <= leafSize) {
      left = NULL; right = NULL;
    } else {
      intT k = findWidest();
      floatT xM = (pMax[k]+pMin[k])/2;

      // Split items by xM (serial)
      intT median = splitItemSerial(xM);

      if (median == 0 || median == size()) {median = ceil(size()/2.0);}

      /* if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) { */
      /*   throw std::runtime_error("Error, kdNode overwrite."); */
      /* } */

      // Recursive construction
      space[0] = nodeT(items.cut(0, median), median, space+1, leafSize);
      space[2*median-1] = nodeT(items.cut(median, size()), size()-median, space+2*median, leafSize);
      left = space;
      right = space+2*median-1;
      left->sib = right;
      right->sib = left;
    }
  }

  template<int _dim, class _objT>
  void kdNode<_dim, _objT>::constructParallel(nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize) {
    boundingBoxParallel();

    sib = NULL;
    if (size() <= leafSize) {
      left = NULL; right = NULL;
    } else {
      intT k = findWidest();
      floatT xM = (pMax[k]+pMin[k])/2;

      // Split items by xM in dim k (parallel)
      parlay::parallel_for(0, size(),
			   [&](intT i) {
			     if (items[i]->at(k)<xM) flags[i]=1;
			     else flags[i] = 0;});
      auto mySplit = kdTreeInternal::split_two(items, flags);
      auto splited = mySplit.first;
      intT median = mySplit.second;
      parlay::parallel_for(0, size(), [&](intT i) {items[i] = splited[i];}); // Copy back

      if (median == 0 || median == size()) {median = (size()/2.0);}

      /* if (!space[0].isEmpty() || !space[2*median-1].isEmpty()) { */
      /*   throw std::runtime_error("Error, kdNode overwrite."); */
      /* } */

      // Recursive construction
      parlay::par_do([&](){space[0] = nodeT(items.cut(0, median), median, space+1, flags.cut(0, median), leafSize);},
		     [&](){space[2*median-1] = nodeT(items.cut(median, size()), size()-median, space+2*median, flags.cut(median, size()), leafSize);});
      left = space;
      right = space+2*median-1;
      left->sib = right;
      right->sib = left;
    }
  }

  template<int _dim, class _objT>
  kdNode<_dim, _objT>::kdNode(parlay::slice<_objT**, _objT**> itemss,
			      intT nn,
			      nodeT *space,
			      parlay::slice<bool*, bool*> flags,
			      intT leafSize): items(itemss)
  {
    resetId();
    if (size()>2000) constructParallel(space, flags, leafSize);
    else constructSerial(space, leafSize);
  }

  template<int _dim, class _objT>
  kdNode<_dim, _objT>::kdNode(parlay::slice<_objT**, _objT**> itemss,
			      intT nn,
			      nodeT *space,
			      intT leafSize): items(itemss)
  {
    resetId();
    constructSerial(space, leafSize);
  }

  template <typename nodeT>
  double nodeDistance(nodeT* n1, nodeT* n2) {
    using floatT = typename nodeT::objT::floatT;

    for (int d = 0; d < n1->dim; ++ d) {
      if (n1->getMin(d) > n2->getMax(d) || n2->getMin(d) > n1->getMax(d)) {
	// disjoint at dim d, and intersect on dim < d
	floatT rsqr = 0;
	for (int dd = d; dd < n1->dim; ++ dd) {
	  floatT tmp = std::max(n1->getMin(dd) - n2->getMax(dd),
				n2->getMin(dd) - n1->getMax(dd));
	  tmp = std::max(tmp, (floatT)0);
	  rsqr += tmp*tmp;
	}
	return sqrt(rsqr);
      }
    }
    return 0; // could be intersecting
  }

  template <typename nodeT>
  double nodeFarDistance(nodeT* n1, nodeT* n2) {
    using floatT = typename nodeT::objT::floatT;
    floatT result = 0;
    for (int d = 0; d < n1->dim; ++ d) {
      floatT tmp = std::max(n1->getMax(d), n2->getMax(d)) - std::min(n1->getMin(d), n2->getMin(d));
      result += tmp * tmp;
    }
    return sqrt(result);
  }

  template<int dim, class objT>
  kdNode<dim, objT>* buildKdTree(parlay::slice<objT*, objT*> P,
				 bool parallel,
				 size_t leafSize,
				 parlay::sequence<objT*>* items) {
    typedef kdNode<dim, objT> nodeT;

    size_t n = P.size();

    bool freeItems = false;
    if (!items) {
      items = new parlay::sequence<objT*>(n);
      freeItems = true;
    }

    parlay::parallel_for(0, n, [&](size_t i) {items->at(i)=&P[i];});

    parlay::slice<objT**, objT**> itemSlice = items->cut(0, items->size());

    auto root = (nodeT*) malloc(sizeof(nodeT)*(2*n-1));

    if (parallel) {
      auto flags = parlay::sequence<bool>(n);
      auto flagSlice = parlay::slice(flags.begin(), flags.end());
      root[0] = nodeT(itemSlice, n, root+1, flagSlice, leafSize);
    } else {
      root[0] = nodeT(itemSlice, n, root+1, leafSize);
    }

    //if (freeItems) free(items);

    return root;
  }

  template<int dim, class objT>
  kdNode<dim, objT>* buildKdTree(parlay::sequence<objT>& P,
				 bool parallel,
				 size_t leafSize,
				 parlay::sequence<objT*>* items) {
    return buildKdTree<dim, objT>(parlay::make_slice(P), parallel, leafSize, items);
  }

} // End namespace pargeo
