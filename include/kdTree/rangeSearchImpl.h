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
#include "parlay/sequence.h"
#include "kdTree.h"
#include "pargeo/point.h"

namespace pargeo::kdTree
{

  template <int dim, typename nodeT, typename objT>
  void rangeHelper(nodeT *tree, objT &q, point<dim> qMin, point<dim> qMax,
                   double radius, parlay::sequence<objT *> &out)
  {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if (relation == tree->boxExclude)
    {
      return;
    }
    else if (relation == tree->boxInclude)
    {
      for (size_t i = 0; i < tree->size(); ++i)
      {
        objT *p = tree->getItem(i);
        if (p->dist(q) <= radius)
          out.push_back(p);
      }
    }
    else
    { // intersect
      if (tree->isLeaf())
      {
        for (size_t i = 0; i < tree->size(); ++i)
        {
          objT *p = tree->getItem(i);
          double dist = q.dist(*p);
          if (dist <= radius)
            out.push_back(p);
        }
      }
      else
      {
        rangeHelper<dim, nodeT, objT>(tree->L(), q, qMin, qMax, radius, out);
        rangeHelper<dim, nodeT, objT>(tree->R(), q, qMin, qMax, radius, out);
      }
    }
  }

  template <int dim, typename objT>
  parlay::sequence<objT *> rangeSearch(
      node<dim, objT> *tree,
      objT query,
      double radius)
  {
    auto out = parlay::sequence<objT *>();
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - radius;
      qMin[i] = tmp;
      qMax[i] = tmp + radius * 2;
    }
    rangeHelper<dim, node<dim, objT>, objT>(tree, query, qMin, qMax,
                                            radius, out);
    return out;
  }

  template <int dim, typename objT>
  parlay::sequence<size_t> bruteforceRange(parlay::sequence<objT> &elems,
                                           objT query,
                                           double radius)
  {
    auto out = parlay::sequence<objT>();
    auto flag = parlay::sequence<size_t>(elems.size(), elems.size());
    parallel_for(0, elems.size(), [&](size_t i)
                 {
                   if (elems[i].dist(query) <= radius)
                     flag[i] = i;
                 });
    return parlay::filter(make_slice(flag), [&](size_t i)
                          { return i < elems.size(); });
  }

  template <int dim, typename nodeT, typename objT>
  void orthRangeHelper(nodeT *tree, point<dim> qMin, point<dim> qMax,
                       parlay::sequence<objT *> &out)
  {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if (relation == tree->boxExclude)
    {
      return;
    }
    else if (relation == tree->boxInclude)
    {
      for (size_t i = 0; i < tree->size(); ++i)
      {
        objT *p = tree->getItem(i);
        out.push_back(p);
      }
    }
    else
    { // intersect
      if (tree->isLeaf())
      {
        for (size_t i = 0; i < tree->size(); ++i)
        {
          objT *p = tree->getItem(i);
          objT _p = *p;
          bool in = true;
          for (int d = 0; d < dim; ++d)
          {
            if (_p[d] > qMax[d] || _p[d] < qMin[d])
              in = false;
          }
          if (in)
            out.push_back(p);
        }
      }
      else
      {
        orthRangeHelper<dim, nodeT, objT>(tree->L(), qMin, qMax, out);
        orthRangeHelper<dim, nodeT, objT>(tree->R(), qMin, qMax, out);
      }
    }
  }

  template <int dim, typename objT>
  parlay::sequence<objT *> orthogonalRangeSearch(node<dim, objT> *tree,
                                                 objT query,
                                                 double halfLen)
  {
    auto out = parlay::sequence<objT *>();
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    orthRangeHelper<dim, node<dim, objT>, objT>(tree, qMin, qMax, out);
    return out;
  }

  template <int dim, typename objT>
  parlay::sequence<size_t> bruteforceOrthRange(parlay::sequence<objT> &A,
                                               objT query,
                                               double halfLen)
  {
    auto out = parlay::sequence<size_t>();
    point<dim> qMin, qMax;
    for (size_t i = 0; i < dim; i++)
    {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    parallel_for(0, A.size(),
                 [&](size_t i)
                 {
                   bool in = true;
                   for (int d = 0; d < dim; ++d)
                   {
                     if (A[i][d] > qMax[d] || A[i][d] < qMin[d])
                       in = false;
                   }
                   if (in)
                     out.push_back(i);
                 });
    return out;
  }

} // End namespace
