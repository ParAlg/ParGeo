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
#include "pargeo/point.h"
#include <tuple>

namespace pargeo::kdTree
{



  /* Kd-tree node */

  template <int _dim, class _objT>
  class node;

  /* Kd-tree construction and destruction */

  template <int dim, class objT>
  node<dim, objT> *build(parlay::slice<objT *, objT *> P,
                         bool parallel = true,
                         size_t leafSize = 16);

  template <int dim, class objT>
  node<dim, objT> *build(parlay::sequence<objT> &P,
                         bool parallel = true,
                         size_t leafSize = 16);

  template <int dim, class objT>
  void del(node<dim, objT> *tree);

  /* Kd-tree knn search */

  template <int dim, class objT>
  parlay::sequence<size_t> batchKnn(parlay::sequence<objT> &queries,
                                    size_t k,
                                    node<dim, objT> *tree = nullptr,
                                    bool sorted = false);

  /* Kd-tree range search */

  template <int dim, typename objT>
  parlay::sequence<objT *> rangeSearch(
      node<dim, objT> *tree,
      objT query,
      double radius);

  template <int dim, typename objT>
  parlay::sequence<objT *> orthogonalRangeSearch(
      node<dim, objT> *tree,
      objT query,
      double halfLen);

  /* Bichromatic closest pair */

  template <typename nodeT>
  std::tuple<typename nodeT::objT *,
             typename nodeT::objT *,
             typename nodeT::objT::floatT>
  bichromaticClosestPair(nodeT *n1, nodeT *n2);

  /* Well-separated pair decomposition */

  template <typename nodeT>
  struct wsp
  {
    nodeT *u;
    nodeT *v;
    wsp(nodeT *uu, nodeT *vv) : u(uu), v(vv) {}
  };

  template <int dim>
  parlay::sequence<wsp<node<dim, point<dim>>>>
  wellSeparatedPairDecomp(node<dim, point<dim>> *tree, double s = 2);



  /********* Implementations *********/

  template <int _dim, class _objT>
  class tree : public node<_dim, _objT>
  {

  private:
    using baseT = node<_dim, _objT>;

    parlay::sequence<_objT *> *allItems;
    node<_dim, _objT> *space;

  public:
    tree(parlay::slice<_objT *, _objT *> _items,
         typename baseT::intT leafSize = 16)
    {

      typedef tree<_dim, _objT> treeT;
      typedef node<_dim, _objT> nodeT;

      // allocate space for children
      space = (nodeT *)malloc(sizeof(nodeT) * (2 * _items.size() - 1));

      // allocate space for a copy of the items
      allItems = new parlay::sequence<_objT *>(_items.size());

      for (size_t i = 0; i < _items.size(); ++i)
        allItems->at(i) = &_items[i];

      // construct self

      baseT::items = allItems->cut(0, allItems->size());

      baseT::resetId();
      baseT::constructSerial(space, leafSize);
    }

    tree(parlay::slice<_objT *, _objT *> _items,
         parlay::slice<bool *, bool *> flags,
         typename baseT::intT leafSize = 16)
    {

      typedef tree<_dim, _objT> treeT;
      typedef node<_dim, _objT> nodeT;

      // allocate space for children
      space = (nodeT *)malloc(sizeof(nodeT) * (2 * _items.size() - 1));

      // allocate space for a copy of the items
      allItems = new parlay::sequence<_objT *>(_items.size());

      parlay::parallel_for(0, _items.size(), [&](size_t i)
                           { allItems->at(i) = &_items[i]; });

      // construct self

      baseT::items = allItems->cut(0, allItems->size());

      baseT::resetId();
      if (baseT::size() > 2000)
        baseT::constructParallel(space, flags, leafSize);
      else
        baseT::constructSerial(space, leafSize);
    }

    ~tree()
    {
      free(space);
      delete allItems;
    }
  };

  template <int _dim, class _objT>
  class node
  {

  protected:
    using intT = int;

    using floatT = double;

    using pointT = _objT;

    using nodeT = node<_dim, _objT>;

    intT id;

    int k;

    pointT pMin, pMax;

    nodeT *left;

    nodeT *right;

    nodeT *sib;

    parlay::slice<_objT **, _objT **> items;

    inline void minCoords(pointT &_pMin, pointT &p)
    {
      for (int i = 0; i < _pMin.dim; ++i)
        _pMin[i] = std::min(_pMin[i], p[i]);
    }

    inline void maxCoords(pointT &_pMax, pointT &p)
    {
      for (int i = 0; i < _pMax.dim; ++i)
        _pMax[i] = std::max(_pMax[i], p[i]);
    }

    void boundingBoxSerial();

    void boundingBoxParallel();

    intT splitItemSerial(floatT xM);

    inline bool itemInBox(pointT pMin1, pointT pMax1, _objT *item)
    {
      for (int i = 0; i < _dim; ++i)
      {
        if (pMax1[i] < item->at(i) || pMin1[i] > item->at(i))
          return false;
      }
      return true;
    }

    inline intT findWidest()
    {
      floatT xM = -1;
      for (int kk = 0; kk < _dim; ++kk)
      {
        if (pMax[kk] - pMin[kk] > xM)
        {
          xM = pMax[kk] - pMin[kk];
          k = kk;
        }
      }
      return k;
    }

    void constructSerial(nodeT *space, intT leafSize);

    void constructParallel(nodeT *space, parlay::slice<bool *, bool *> flags, intT leafSize);

  public:
    using objT = _objT;

    static constexpr int dim = _dim;

    inline nodeT *L() { return left; }

    inline nodeT *R() { return right; }

    inline nodeT *siblin() { return sib; }

    inline intT size() { return items.size(); }

    inline _objT *operator[](intT i) { return items[i]; }

    inline _objT *at(intT i) { return items[i]; }

    inline bool isLeaf() { return !left; }

    inline _objT *getItem(intT i) { return items[i]; }

    inline pointT getMax() { return pMax; }

    inline pointT getMin() { return pMin; }

    inline floatT getMax(int i) { return pMax[i]; }

    inline floatT getMin(int i) { return pMin[i]; }

    // inline void setEmpty() { id = -2; }

    // inline bool isEmpty() { return id == -2; }

    inline bool hasId()
    {
      return id != -1;
    }

    inline void setId(intT _id)
    {
      id = _id;
    }

    inline void resetId()
    {
      id = -1;
    }

    inline intT getId()
    {
      return id;
    }

    static const int boxInclude = 0;

    static const int boxOverlap = 1;

    static const int boxExclude = 2;

    inline floatT diag()
    {
      floatT result = 0;
      for (int d = 0; d < _dim; ++d)
      {
        floatT tmp = pMax[d] - pMin[d];
        result += tmp * tmp;
      }
      return sqrt(result);
    }

    inline floatT lMax()
    {
      floatT myMax = 0;
      for (int d = 0; d < _dim; ++d)
      {
        floatT thisMax = pMax[d] - pMin[d];
        if (thisMax > myMax)
        {
          myMax = thisMax;
        }
      }
      return myMax;
    }

    inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2)
    {
      bool exclude = false;
      bool include = true; //1 include 2
      for (int i = 0; i < _dim; ++i)
      {
        if (pMax1[i] < pMin2[i] || pMin1[i] > pMax2[i])
          exclude = true;
        if (pMax1[i] < pMax2[i] || pMin1[i] > pMin2[i])
          include = false;
      }
      if (exclude)
        return boxExclude;
      else if (include)
        return boxInclude;
      else
        return boxOverlap;
    }

    node();

    node(parlay::slice<_objT **, _objT **> itemss,
         intT nn,
         nodeT *space,
         parlay::slice<bool *, bool *> flags,
         intT leafSize = 16);

    node(parlay::slice<_objT **, _objT **> itemss,
         intT nn,
         nodeT *space,
         intT leafSize = 16);

    virtual ~node()
    {
    }
  };

} // End namespace pargeo::kdTree

#include "treeImpl.h"
#include "knnImpl.h"
#include "rangeSearchImpl.h"
#include "bccpImpl.h"
#include "wspdImpl.h"
