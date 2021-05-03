#pragma once

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "kdTree.h"
#include "point.h"

namespace pargeo {
  using namespace parlay;

  template<int dim, typename nodeT, typename objT>
    void rangeHelper(nodeT* tree, objT& q, point<dim> qMin, point<dim> qMax,
			double radius, sequence<size_t>& out, objT* A) {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if(relation == tree->boxExclude) {
      return;
    } else if (relation == tree->boxInclude) {
      for (size_t i = 0; i < tree->size(); ++i) {
	objT* p = tree->getItem(i);
	out.push_back(p - A);
      }
    } else { // intersect
      if (tree->isLeaf()) {
	for (size_t i = 0; i < tree->size(); ++ i) {
	  objT* p = tree->getItem(i);
	  double dist = q.dist(*p);
	  if (dist <= radius)
	    out.push_back(p - A);
	}
      } else {
	rangeHelper<dim, nodeT, objT>(tree->L(), q, qMin, qMax, radius, out, A);
	rangeHelper<dim, nodeT, objT>(tree->R(), q, qMin, qMax, radius, out, A);
      }
    }
  }

  template<int dim, typename objT>
  sequence<size_t> kdTreeRange(sequence<objT>& A, kdNode<dim, objT>* tree,
			       objT query, double radius) {
    auto out = parlay::sequence<size_t>();
    point<dim> qMin, qMax;
    for (size_t i=0; i<dim; i++) {
      auto tmp = query[i] - radius;
      qMin[i] = tmp;
      qMax[i] = tmp + radius * 2;
    }
    rangeHelper<dim, kdNode<dim, objT>, objT>(tree, query, qMin, qMax,
						 radius, out, A.data());
    return out;
  }

  template<int dim, typename objT>
  sequence<size_t> bruteforceRange(sequence<objT>& elems, objT query, double radius) {
    auto out = parlay::sequence<objT>();
    auto flag = parlay::sequence<size_t>(elems.size(), elems.size());
    parallel_for(0, elems.size(), [&](size_t i) {
	if (elems[i].dist(query) <= radius)
	  flag[i] =  i;
      });
    return parlay::filter(make_slice(flag), [&](size_t i) {
	return i < elems.size();
      });
  }

  template<int dim, typename nodeT, typename objT>
    void orthRangeHelper(nodeT* tree, point<dim> qMin, point<dim> qMax,
			 sequence<size_t>& out, objT* A) {
    int relation = tree->boxCompare(qMin, qMax, tree->getMin(), tree->getMax());

    if(relation == tree->boxExclude) {
      return;
    } else if (relation == tree->boxInclude) {
      for (size_t i = 0; i < tree->size(); ++i) {
	objT* p = tree->getItem(i);
	out.push_back(p - A);
      }
    } else { // intersect
      if (tree->isLeaf()) {
	for (size_t i = 0; i < tree->size(); ++ i) {
	  objT *p = tree->getItem(i);
	  objT _p = *p;
	  bool in = true;
	  for (int d = 0; d < dim; ++ d) {
	    if (_p[d] > qMax[d] || _p[d] < qMin[d])
	      in = false;
	  }
	  if (in)
	    out.push_back(p - A);
	}
      } else {
	orthRangeHelper<dim, nodeT, objT>(tree->L(), qMin, qMax, out, A);
	orthRangeHelper<dim, nodeT, objT>(tree->R(), qMin, qMax, out, A);
      }
    }
  }

  template<int dim, typename objT>
  sequence<size_t> kdTreeOrthRange(sequence<objT>& A, kdNode<dim, objT>* tree,
			       objT query, double halfLen) {
    auto out = parlay::sequence<size_t>();
    point<dim> qMin, qMax;
    for (size_t i=0; i<dim; i++) {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    orthRangeHelper<dim, kdNode<dim, objT>, objT>(tree, qMin, qMax,
						 out, A.data());
    return out;
  }

  template<int dim, typename objT>
  sequence<size_t> bruteforceOrthRange(sequence<objT>& A,
				       objT query, double halfLen) {
    auto out = parlay::sequence<size_t>();
    point<dim> qMin, qMax;
    for (size_t i=0; i<dim; i++) {
      auto tmp = query[i] - halfLen;
      qMin[i] = tmp;
      qMax[i] = tmp + halfLen * 2;
    }
    parallel_for(0, A.size(),
		 [&](size_t i) {
		   bool in = true;
		   for (int d = 0; d < dim; ++ d) {
		     if (A[i][d] > qMax[d] || A[i][d] < qMin[d])
		       in = false;
		   }
		   if (in) out.push_back(i);
		 });
    return out;
  }

} // End namespace
