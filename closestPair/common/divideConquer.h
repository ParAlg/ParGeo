// This code is part of the project "A Parallel Batch-Dynamic Data Structure
// for the Closest Pair Problem"
// Copyright (c) 2020 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
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

#ifndef DIVIDE_CONQUER_H
#define DIVIDE_CONQUER_H

#include "geometry.h"
#include "shared.h"
#include "pbbs/parallel.h"
#include "pbbs/sampleSort.h"
#include "pbbs/quickSort.h"
#include "pbbs/sequence.h"
using namespace std;;

// *************************************************************
//    Utils
// *************************************************************

//special comparator that compares a point along a certain dimension
template<int dim>
bool comparePoint(const point<dim> i, const point<dim> j, int k) {
  return i.x[k] < j.x[k];
}

template<int dim>
struct findSlab {
  double xm, delta;
  int k; //dimension along which the slab is computed
  findSlab(double xmm, double deltaa, int kk): xm(xmm), delta(deltaa), k(kk) {
    if (k >= dim) {
      cout << "findSlab exceeds point dimension, abort" << endl; abort(); }
  }
  bool operator() (const point<dim> i) {
    if (i.x[k] <= xm+delta && i.x[k] >= xm-delta) return true;
    else return false;
  }
};

inline intT intPow(intT base, intT exp) {
  intT r = 1;
  for (intT i=0; i<exp; ++i) r *= base;
  return r;
}

template<int dim>
inline intT numToCheck() {return intPow(2, dim*2-1)-1;}

// *************************************************************
//    Divide conquer
// *************************************************************

/**
 * Sequentially, computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point array.
 * @param n length of P.
 * @param Pp pre-allocated auxiliary memory of length n.
 * @param k dimensionality of the call, k = [0,dim)
 * @return the closest pair
 */
template<int dim>
pointPair<dim> divideConquerSerial(point<dim>* P, intT n, point<dim>* Pp, int k) {
  if (n < 100) return bruteForceSerial<dim>(P, n);

  auto pLess = [&] (point<dim> a, point<dim> b) {
                 return comparePoint<dim>(a, b, k);};
  quickSortSerial(P, n, pLess);
  intT xm = intT(n/2);
  double xmm = P[xm].x[k];

  point<dim> *L = P;
  point<dim> *R = P+xm;
  point<dim> *Lp = Pp;
  point<dim> *Rp = Pp+xm;
  pointPair<dim> Lr = divideConquerSerial<dim>(L, xm, Lp, k);
  pointPair<dim> Rr = divideConquerSerial<dim>(R, n-xm, Rp, k);
  pointPair<dim> better = Lr.dist < Rr.dist ? Lr : Rr;
  if (better.dist == 0) return better;

  //create slab in dimension k, of size np and stored in aux-memory Pp
  intT np = sequence::filterSerial(P, Pp, n, findSlab<dim>(xmm, better.dist, k));

  pointPair<dim> slabR = pointPair<dim>();

  if (k >= 2) {
    //recursive case, solve the slab in one dimension lower
    if (np > n/2-1) {
      //does not have sufficient auxiliary memory
      point<dim>* Ppp = newA(point<dim>, np);
      slabR = divideConquerSerial<dim>(Pp, np, Ppp, k-1);
      free(Ppp);
    } else {
      slabR = divideConquerSerial<dim>(Pp, np, Pp+np, k-1);
    }
  } else if (k == 1) {
    // base case
    auto pLess = [&] (point<dim> a, point<dim> b)
                 {return comparePoint<dim>(a, b, 0);};
    quickSortSerial(Pp, np, pLess);
    const intT check = numToCheck<dim>();
    for(intT i = 0; i < np; ++ i) { //for each point in the slab
      for (intT j = 0; j < check; ++ j) { //check constant number of points after
        if (i+j+1<np) {
          slabR.closer(Pp[i], Pp[i+j+1]);}
      }}
  } else {
    //should never reach here
    cout << "wrong k = " << k << ", abort" << endl; abort();
  }
  if (slabR.dist < better.dist) return slabR;
  else return better;
}

/**
 * In parallel, computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point array.
 * @param n length of P.
 * @param Pp pre-allocated auxiliary memory of length n.
 * @param k dimensionality of the call, k = [0,dim)
 * @return the closest pair
 */
template<int dim>
pointPair<dim> divideConquerParallel(point<dim>* P, intT n, point<dim>* Pp, int k) {
  if (n < 10000) return divideConquerSerial(P, n, Pp, k);

  auto pLess = [&] (point<dim> a, point<dim> b)
               {return comparePoint<dim>(a, b, k);};
  sampleSort(P, n, pLess);
  intT xm = intT(n/2);
  double xmm = P[xm].x[k];

  point<dim> *L = P;
  point<dim> *R = P+xm;
  point<dim> *Lp = Pp;
  point<dim> *Rp = Pp+xm;
  pointPair<dim> Lr = cilk_spawn divideConquerParallel<dim>(L, xm, Lp, k);
  pointPair<dim> Rr = divideConquerParallel<dim>(R, n-xm, Rp, k);
  cilk_sync;

  pointPair<dim> better = Lr.dist < Rr.dist ? Lr : Rr;
  if (better.dist == 0) return better;

  //create slab in dimension k, of size np and stored in aux-memory Pp
  intT np = sequence::filter(P, Pp, n, findSlab<dim>(xmm, better.dist, k));

  pointPair<dim> slabR = pointPair<dim>(P[0], P[1], intMax());;
  if (k >= 2) {
    //recursive case, solve the slab in one dimension lower
    if (np > n/2-1) {
      //does not have sufficient auxiliary memory
      point<dim>* Ppp = newA(point<dim>, np);
      slabR = divideConquerParallel<dim>(Pp, np, Ppp, k-1);
      free(Ppp);
    } else {
      slabR = divideConquerParallel<dim>(Pp, np, Pp+np, k-1);
    }
  } else if (k == 1) {
    //base case
    auto pLess = [&] (point<dim> a, point<dim> b) {
                   return comparePoint<dim>(a, b, 0);};
    const intT check = numToCheck<dim>();
    if (np > 100) {
      //go parallel
      sampleSort(Pp, np, pLess);
      pointPair<dim>* A = newA(pointPair<dim>, np*check);
      par_for(intT i=0; i<np; ++i) { //for each point in the slab in parallel
        for (intT j = 0; j < check; ++ j) { //check constant number of points after
          if (i+j+1<np) {
            A[i*check+j] = pointPair<dim>(Pp[i], Pp[i+j+1]);
          } else {
            A[i*check+j] = pointPair<dim>(Pp[i], Pp[i+j+1], intMax());
          }
        }
      }
      auto getDist = [&](intT i) {return A[i].dist;};
      intT I = sequence::minIndex<double, intT>(0, np*check, getDist);
      slabR = A[I];
    } else {
      //go serial
      quickSortSerial(Pp, np, pLess);
      for(intT i = 0; i < np; ++ i) { //for each point in the slab
        for (intT j = 0; j < check; ++ j) { //check constant number of points after
          if (i+j+1<np) {
            slabR.closer(Pp[i], Pp[i+j+1]);}
        }}
    }
  } else {
    cout << "wrong k = " << k << ", abort" << endl; abort();
  }
  if (slabR.dist < better.dist) return slabR;
  else return better;
}

/**
 * Computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point array.
 * @param n length of P.
 * @param serial optional boolean variable which should be set to true if wish to run serialy.
 * @return the closest pair
 */
template<int dim>
pointPair<dim> divideConquer(point<dim>* P, intT n, bool serial = false) {
  point<dim>* Pp = newA(point<dim>, n*2);
  par_for(intT i=0; i<n; ++i) {Pp[i]=P[i];}
  pointPair<dim> r;
  if (serial) r = divideConquerSerial<dim>(Pp, n, Pp+n, dim-1);
  else r = divideConquerParallel<dim>(Pp, n, Pp+n, dim-1);
  free(Pp);
  return r;
}

#endif
