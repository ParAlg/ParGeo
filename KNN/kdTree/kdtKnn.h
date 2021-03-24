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

#ifndef KDT_KNN_H
#define KDT_KNN_H

#include <limits> // numeric_limits
#include <algorithm> // nth_element
#include "parlay/parallel.h"
#include "parlay/sequence.h"

namespace knnBuf {

  typedef int intT;
  typedef double floatT;

  template <typename T>
  struct elem {
    floatT cost;// Non-negative
    T entry;
    elem(floatT t_cost, T t_entry) : cost(t_cost), entry(t_entry) {}
    elem() : cost(std::numeric_limits<floatT>::max()) {}
    bool operator<(const elem& b) const {
      if (cost < b.cost) return true;
      return false;}
  };

  template <typename T>
  struct buffer {
    typedef parlay::slice<elem<T>*, elem<T>*> sliceT;
    intT k;
    intT ptr;
    sliceT buf;

    buffer(intT t_k, sliceT t_buf): k(t_k), ptr(0), buf(t_buf) {}

    inline void reset() {ptr = 0;}

    bool hasK() {return ptr >= k;}

    elem<T> keepK() {
      if (ptr < k) abort();
      ptr = k;
      std::nth_element(buf.begin(), buf.begin()+k, buf.end());
      return buf[k];
    }

    void insert(elem<T> t_elem) {
      buf[ptr++] = t_elem;
      if (ptr >= buf.size()) keepK();
    }

    elem<T> operator[](intT i) {
      if (i < ptr) return buf[i];
      else return elem<T>();
    }
  };
}

namespace kdtKnn {

  using namespace knnBuf;

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
	cout << "error, kdTree splitting singleton, abort" << endl;abort();}
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

    inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2) {
      bool exclude = false;
      bool include = true;//1 include 2
      for(int i=0; i<dim; ++i) {
	if (pMax1[i]<pMin2[i] || pMin1[i]>pMax2[i]) exclude = true;
	if (pMax1[i]<pMax2[i] || pMin1[i]>pMin2[i]) include = false;
      }
      if (exclude) return boxExclude;
      else if (include) return boxInclude;
      else return boxOverlap;
    }

    inline bool itemInBox(pointT pMin1, pointT pMax1, objT* item) {
      for(int i=0; i<dim; ++i) {
	if (pMax1[i]<item->at(i) || pMin1[i]>item->at(i)) return false;
      }
      return true;
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
	  cout << "error, kdNode overwrite, abort" << endl;abort();}

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
	  cout << "error, kdNode overwrite, abort" << endl;abort();}

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

    void knnRangeHelper(objT& q, pointT qMin, pointT qMax, floatT radius, buffer<objT*>& out);
    void knnRange(objT& q, floatT radius, buffer<objT*>& out);
    void knnHelper(objT& q, buffer<objT*>& out);
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

  template<int dim, class objT>
  void kdNode<dim, objT>::knnRangeHelper(objT& q, pointT qMin, pointT qMax, floatT radius, buffer<objT*>& out) {
    int relation = boxCompare(qMin, qMax, getMin(), getMax());

    if(relation == boxExclude) {
      return;
    } else if (relation == boxInclude) {
      for (intT i=0; i<size(); ++i) {
	objT* p = getItem(i);
	out.insert(elem(q.dist(*p), p));
      }
    } else { // intersect
      if (isLeaf()) {
	for (intT i=0; i < size(); ++ i) {
	  objT* p = getItem(i);
	  float dist = q.dist(*p);
	  if (dist <= radius) {out.insert(elem(dist, p));}
	}
      } else {
	L()->kdNode<dim, objT>::knnRangeHelper(q, qMin, qMax, radius, out);
	R()->kdNode<dim, objT>::knnRangeHelper(q, qMin, qMax, radius, out);
      }
    }
  }

  template<int dim, class objT>
  void kdNode<dim, objT>::knnRange(objT& q, floatT radius, buffer<objT*>& out) {
    pointT qMin, qMax;
    for (intT i=0; i<dim; i++) {
      auto tmp = q[i]-radius;
      qMin[i] = tmp;
      qMax[i] = tmp+radius*2;
    }
    kdNode<dim, objT>::knnRangeHelper(q, qMin, qMax, radius, out);
  }

  template<int dim, class objT>
  void kdNode<dim, objT>::knnHelper(objT& q, buffer<objT*>& out) {
    // find the leaf first
    int relation = boxCompare(getMin(), getMax(), pointT(q.coords()), pointT(q.coords()));
    if (relation == boxExclude) {
      return;
    } else {
      if (isLeaf()) {
	// basecase
	for (intT i=0; i<size(); ++ i) {
	  objT* p = getItem(i);
	  out.insert(elem(q.dist(*p), p));}
      } else {
	L()->kdNode<dim, objT>::knnHelper(q, out);
	R()->kdNode<dim, objT>::knnHelper(q, out);
      }
    }

    if (!out.hasK()) {
      if (siblin() == NULL) {
	cout << "error, knnHelper reached root node without enough neighbors. k = " << k << endl;
	abort();
      }
      for (intT i=0; i<siblin()->size(); ++i) {
	objT* p = siblin()->getItem(i);
	out.insert(elem(q.dist(*p), p));}
    } else { // Buffer filled to a least k
      if (siblin() != NULL) {
	elem tmp = out.keepK();
	siblin()->kdNode<dim, objT>::knnRange(q, tmp.cost, out);}
    }
  }

  template<int dim, class objT>
  parlay::sequence<size_t> kdtKnn(parlay::sequence<pargeo::point<dim>> &queries, size_t k) {
    kdNode<dim, pargeo::point<dim>>* tree = buildKdt<dim, pargeo::point<dim>>(queries, true);
    auto out = parlay::sequence<elem<pargeo::point<dim>*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());
    parlay::parallel_for(0, queries.size(), [&](intT i) {
					      buffer buf = buffer<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      tree->knnHelper(queries[i], buf);
					      buf.keepK();
					      for(intT j=0; j<k; ++j) {
						idx[i*k+j] = buf[j].entry - queries.data();
						//cout << buf[j].cost << endl;
					      }
					      //cout << endl;
					    });
    free(tree);
    return idx;
  }

  template<int dim, class objT>
  parlay::sequence<size_t> bruteforceKnn(parlay::sequence<pargeo::point<dim>> &queries, size_t k) {
    auto out = parlay::sequence<elem<pargeo::point<dim>*>>(2*k*queries.size());
    auto idx = parlay::sequence<size_t>(k*queries.size());
    parlay::parallel_for(0, queries.size(), [&](intT i) {
					      objT q = queries[i];
					      buffer buf = buffer<objT*>(k, out.cut(i*2*k, (i+1)*2*k));
					      for(intT j=0; j<queries.size(); ++j) {
						objT* p = &queries[j];
						buf.insert(elem(q.dist(p), p));
					      }
					      buf.keepK();
					      for(intT j=0; j<k; ++j) {
						idx[i*k+j] = buf[j].entry - queries.data();
						//cout << buf[j].cost << endl;
					      }
					      //cout << endl;
					    });
    return idx;
  }
}

#endif
