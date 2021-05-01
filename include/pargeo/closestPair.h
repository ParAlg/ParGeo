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

#pragma once

#include <tuple>
#include <limits>
#include "pargeo/point.h"
#include "parlay/parallel.h"
#include "parlay/internal/sample_sort.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
static constexpr double floatMax = std::numeric_limits<double>::max();

template<int dim>
struct pointPair {
  using pt = point<dim>;
  using pp = pointPair<dim>;
  using floatT = typename pt::floatT;

  floatT dist;

  pt u,v;

  pointPair(): dist(floatMax) {u=pt();v=pt();}

  pointPair(pt uu, pt vv, floatT distt): u(uu), v(vv), dist(distt) {}

  pointPair(pt uu, pt vv): u(uu), v(vv) {
    dist = u.dist(v);}

  bool isEmpty() {return dist==floatMax;}

  void closer(pt uu, pt vv) {
    floatT distt = uu.dist(vv);
    if (distt<dist) {
      u = uu;
      v = vv;
      dist = distt;}
  }

  void closer(pt uu, pt vv, pt distt) {
    if (distt<dist) {
      u = uu;
      v = vv;
      dist = distt;}
  }

  void closer(pt pair) {
    if (pair.dist < dist) {
      u = pair.u;
      v = pair.v;
      dist = pair.dist;}
  }

  friend bool operator<(pp a, pp b) {
    if (a.dist == b.dist) {
      for (int i=0; i<dim; ++i) {
        if(a.u[i]==b.u[i]) continue;
        else return a.u[i]<b.u[i];}
      return false;
    }
    else return a.dist < b.dist;
  }
};

template<int dim>
pointPair<dim> bruteForceParallel(slice<point<dim>*,point<dim>*> P) {
  auto A = sequence<pointPair<dim>>(P.size()*P.size());
  parallel_for(0, P.size(), [&](size_t i) {
      for(size_t j = 0; j < P.size(); ++ j) {
	if (i == j) {
	  A[i+j*P.size()] = pointPair<dim>(P[i], P[j], floatMax);
	} else {
	  A[i+j*P.size()] = pointPair<dim>(P[i], P[j]);
	}
      }
    });
  auto R = *min_element(make_slice(A), [&](pointPair<dim> i, pointPair<dim> j) { return i.dist < j.dist; });
  return R;
}

template<int dim>
pointPair<dim> bruteForceSerial(slice<point<dim>*,point<dim>*> P) {
  pointPair<dim> A = pointPair<dim>(P[0], P[1], floatMax);
  parallel_for(0, P.size(), [&](size_t i) {
      for(size_t j = i+1; j < P.size(); ++ j) {
        A.closer(P[i], P[j]);
      }
  });
  return A;
}

// *************************************************************
//    Helper functions
// *************************************************************

// todo change double to appropriate type
template<int dim>
struct findSlab {
  double xm, delta;
  int k; //dimension along which the slab is computed
  findSlab(double xmm, double deltaa, int kk): xm(xmm), delta(deltaa), k(kk) {
    if (k >= dim) {
      throw std::runtime_error("findSlab exceeds point dimension, abort"); }
  }
  bool operator() (const point<dim> i) {
    if (i.x[k] <= xm+delta && i.x[k] >= xm-delta) return true;
    else return false;
  }
};

inline size_t intPow(size_t base, size_t exp) {
  size_t r = 1;
  for (size_t i=0; i<exp; ++i) r *= base;
  return r;
}

template<int dim>
inline size_t numToCheck() {return intPow(2, dim*2-1)-1;}

// *************************************************************
//    Divide conquer
// *************************************************************

/**
 * Sequentially, computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point sequence.
 * @param k dimensionality of the call, k = [0,dim)
 * @return the closest pair
 */
template<int dim>
pointPair<dim> divideConquerSerial(slice<point<dim>*,point<dim>*> P, int k) {
  using pt = point<dim>;
  using pp = pointPair<dim>;
  if (P.size() < 100) return bruteForceSerial<dim>(P);

  internal::seq_sort_inplace(P, [&](pt i, pt j){return i[k] < j[k];}, false);
  size_t xm = size_t(P.size()/2);
  double xmm = P[xm].x[k];

  auto L = P.cut(0, xm);
  auto R = P.cut(xm, P.size());
  pp Lr = divideConquerSerial<dim>(L, k);
  pp Rr = divideConquerSerial<dim>(R, k);
  pp better = Lr.dist < Rr.dist ? Lr : Rr;
  if (better.dist == 0) return better;

  //create slab in dimension k, of size np and stored in aux-memory Pp
  auto Pp = filter(P, findSlab<dim>(xmm, better.dist, k));// serial todo

  auto slabR = pp();

  if (k >= 2) {
    slabR = divideConquerSerial<dim>(make_slice(Pp), k-1);
  } else if (k == 1) {
    internal::seq_sort_inplace(make_slice(Pp), [&](pt i, pt j){return i[0] < j[0];}, false);
    const size_t check = numToCheck<dim>();
    for(size_t i = 0; i < Pp.size(); ++ i) { //for each point in the slab
      for (size_t j = 0; j < check; ++ j) { //check constant number of points after
        if (i+j+1<Pp.size()) {
          slabR.closer(Pp[i], Pp[i+j+1]);}
      }}
  } else {
    throw std::runtime_error("error, wrong k, abort");
  }
  if (slabR.dist < better.dist) return slabR;
  else return better;
}

/**
 * In parallel, computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point sequence.
 * @param k dimensionality of the call, k = [0,dim)
 * @return the closest pair
 */
template<int dim>
pointPair<dim> divideConquerParallel(slice<point<dim>*,point<dim>*> P, int k) {
  using pt = point<dim>;
  using pp = pointPair<dim>;
  if (P.size() < 10000) return divideConquerSerial(P, k);

  sort_inplace(P, [&](pt i, pt j){return i[k] < j[k];});
  size_t xm = size_t(P.size()/2);
  double xmm = P[xm].x[k];

  auto L = P.cut(0, xm);
  auto R = P.cut(xm, P.size());
  pp Lr, Rr;
  par_do([&](){Lr = divideConquerParallel<dim>(L, k);},
	 [&](){Rr = divideConquerParallel<dim>(R, k);});
  pp better = Lr.dist < Rr.dist ? Lr : Rr;
  if (better.dist == 0) return better;

  //create slab in dimension k, of size np and stored in aux-memory Pp
  auto Pp = filter(P, findSlab<dim>(xmm, better.dist, k));
  size_t np = Pp.size();

  auto slabR = pp(P[0], P[1], floatMax);
  if (k >= 2) {
    //recursive case, solve the slab in one dimension lower
    slabR = divideConquerParallel<dim>(make_slice(Pp), k-1);
  } else if (k == 1) {
    //base case
    const size_t check = numToCheck<dim>();
    if (np > 100) {
      //go parallel
      sort_inplace(make_slice(Pp), [&](pt i, pt j){return i[0] < j[0];});
      auto A = sequence<pp>(np*check);
      parallel_for(0, np, [&](size_t i) { //for each point in the slab in parallel
	  for (size_t j = 0; j < check; ++ j) { //check constant number of points after
	    if (i+j+1<np) {
	      A[i*check+j] = pp(Pp[i], Pp[i+j+1]);
	    } else {
	      A[i*check+j] = pp(Pp[i], Pp[i+j+1], floatMax);
	    }
	  }
	});
      slabR = *min_element(make_slice(A), [&](pp i, pp j) { return i.dist < j.dist; });
    } else {
      //go serial
      internal::seq_sort_inplace(make_slice(Pp), [&](pt i, pt j){return i[0] < j[0];}, false);
      const size_t check = numToCheck<dim>();
      for(size_t i = 0; i < Pp.size(); ++ i) { //for each point in the slab
	for (size_t j = 0; j < check; ++ j) { //check constant number of points after
	  if (i+j+1<Pp.size()) {
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
 * @param serial optional boolean variable which should be set to true if wish to run serialy.
 * @return the closest pair
 */
template<int dim>
pointPair<dim> closestPairDC(sequence<point<dim>>& P, bool serial = false) {
  using pt = point<dim>;
  auto Pp = sequence<pt>(P);
  pointPair<dim> r;

  if (serial) r = divideConquerSerial<dim>(make_slice(Pp), dim-1);
  else r = divideConquerParallel<dim>(make_slice(Pp), dim-1);

  cout << "dist = " << r.dist << endl;

  //cout << "brute = " << bruteForceSerial<dim>(make_slice(P)).dist << endl;
  //cout << "brute = " << bruteForceParallel<dim>(make_slice(P)).dist << endl;
  return r;
}
