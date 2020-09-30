#ifndef PARALLEL_KRUSKAL_H
#define PARALLEL_KRUSKAL_H

// (adapted version)
// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#include <iostream>
#include <limits.h>
#include "pbbs/sequence.h"
#include "pbbs/graph.h"
#include "pbbs/parallel.h"
#include "pbbs/gettime.h"
#include "pbbs/speculative_for.h"
#include "pbbs/unionFind.h"

using namespace std;

#include "pbbs/sampleSort.h"

namespace parKruskal {
  // the par kruskal here work on any array of edges having (u,v,weight)
  // just define an edge wrapper below

  template<class bcpT>
  struct bcpWrap {
    bcpT *E;
    bcpWrap(bcpT *Ein) : E(Ein) {}
    inline intT GetU(intT i) {return E[i].u;}
    inline intT GetV(intT i) {return E[i].v;}
    inline double GetWeight(intT i) {return E[i].dist;}
  };

  // **************************************************************
  //    PARALLEL VERSION OF KRUSKAL'S ALGORITHM
  // **************************************************************

  struct indexedEdge {
    intT u;
    intT v;
    intT id;
  indexedEdge(intT _u, intT _v, intT _id) : u(_u), v(_v), id(_id) {}
  };

  template <typename UFType>
  struct UnionFindStep {
    intT u;  intT v;
    indexedEdge *E;  reservation *R;  UFType* UF;  bool *inST;
    UnionFindStep(indexedEdge* _E, UFType* _UF, reservation* _R, bool* ist)
      : E(_E), R(_R), UF(_UF), inST(ist) {}

    bool reserve(intT i) {
      u = UF->find(E[i].u);
      v = UF->find(E[i].v);
      if (u != v) {
	R[v].reserve(i);
	R[u].reserve(i);
	return 1;
      } else return 0;
    }

    bool commit(intT i) {
      if (R[v].check(i)) {
        R[UF->link(E[i].v, E[i].u)].checkReset(i);
        inST[E[i].id] = 1;
        return 1;
      } else if (R[u].check(i)) {
        R[UF->link(E[i].u, E[i].v)].checkReset(i);
        inST[E[i].id] = 1;
        return 1;
      }
      else {
	return 0;
      }
    }
  };

  struct edgeLess {
    bool operator() (pair<double, intT> a, pair<double, intT> b) {
      return (a.first == b.first) ? (a.second < b.second)
	: (a.first < b.first);}};

  template <typename EdgeWrapper>
  struct myGet {
    EdgeWrapper *A;
    myGet(EdgeWrapper *A) : A(A) {}
    pair<double, intT> operator() (intT i) {
      return pair<double, intT>(A->GetWeight(i),i);
    }
  };

  template <class EdgeWrapper, class F>
  inline intT almostKth(EdgeWrapper *A, pair<double, intT>* B, intT k, intT n, F f) {
    if (n == k) {
      par_for (intT i=0; i < n; i++) {
	B[i] = pair<double, intT>(A->GetWeight(i),i);
      }
      return n;
    }

    intT ssize = min<intT>(1000,n);
    intT stride = n/ssize;
    intT km = (intT)(k * ((double) ssize) / n);
    auto T = (pair<double, intT>*)malloc(ssize*sizeof(pair<double, intT>));
    for (intT i = 0; i < ssize; i++) {
      T[i] = pair<double, intT>(A->GetWeight(i*stride),i*stride);
    }
    sort(T,T+ssize,f);
    pair<double, intT> p = T[km];
    free(T);
    bool *flags = newA(bool,n);
    par_for (intT i=0; i < n; i++) {
      flags[i] = !f(pair<double, intT>(A->GetWeight(i),i),p);
    }
    intT l = sequence::split(B,flags,0,n,myGet<EdgeWrapper>(A));
    free(flags);
    return l;
  }

  //inline pair<intT*,intT> mst(wghEdgeArray<intT> G) {
  template <class EdgeWrapper, class UFType>
  inline intT mst(EdgeWrapper E, intT n, intT m, UFType *UF) {

    //intT l = min<intT>(5*n/4,m); // sort size
    intT l = 0.5*m; // sort size

    auto y = (pair<double, intT>*)malloc(m*sizeof(pair<double, intT>));
    l = almostKth(&E, y, l, m, edgeLess());

    sampleSort(y, l, edgeLess());

    // initialize size n reservation stations
    reservation *R = newA(reservation,n);
    par_for (size_t i = 0; i < n; i++) {
      new ((void*) (R+i)) reservation;}

    //assign each edge an index
    indexedEdge* z = newA(indexedEdge,m);
    par_for (intT i=0; i < l; i++) {
      intT j = y[i].second;
      z[i] = indexedEdge(E.GetU(j),E.GetV(j),j);
    }
    //nextTime("copy to edges");

    bool *mstFlags = newA(bool, m);
    par_for (intT i=0; i < m; i++) mstFlags[i] = 0;
    UnionFindStep<UFType> UFStep(z, UF, R, mstFlags);
    speculative_for(UFStep, 0, l, 8);
    free(z);

    bool *flags = newA(bool,m-l);
    par_for (intT i = 0; i < m-l; i++) {
      intT j = y[i+l].second;
      intT u = UF->find(E.GetU(j));
      intT v = UF->find(E.GetV(j));
      if (u != v) flags[i] = 1;
      else flags[i] = 0;
    }
    auto x = (pair<double, intT>*) malloc(sizeof(pair<double, intT>) * (m-l));
    // input the edges after the l_th, put into x
    intT k = sequence::pack(y+l, x, flags, m-l);
    free(flags);
    free(y);

    sampleSort(x, k, edgeLess());

    z = newA(indexedEdge, k);
    par_for (intT i=0; i < k; i++) {
      intT j = x[i].second;
      z[i] = indexedEdge(E.GetU(j),E.GetV(j),j);
    }
    free(x);

    UFStep = UnionFindStep<UFType>(z, UF, R, mstFlags);
    speculative_for(UFStep, 0, k, 8);

    free(z);

    intT* mst = newA(intT, m);
    intT nInMst = sequence::packIndex(mst, mstFlags, m);
    free(mstFlags);

    //UF.del();
    free(R);
    return nInMst;
  }
}


#endif
