#pragma once

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "kdTree.h"
#include "point.h"

namespace pargeo {
  using namespace parlay;

  template<int dim, typename nodeT, typename objT>
    void knnRangeHelper(nodeT* tree, objT& q, point<dim> qMin, point<dim> qMax,
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
	knnRangeHelper<dim, nodeT, objT>(tree->L(), q, qMin, qMax, radius, out, A);
	knnRangeHelper<dim, nodeT, objT>(tree->R(), q, qMin, qMax, radius, out, A);
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
    knnRangeHelper<dim, kdNode<dim, objT>, objT>(tree, query, qMin, qMax,
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

}
