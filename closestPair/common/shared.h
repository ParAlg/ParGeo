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

#ifndef SHARED_H
#define SHARED_H

#include "geometry.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"

// *************************************************************
//    Point pair
// *************************************************************

template<int dim>
struct pointPair {
  typedef point<dim> pointT;
  typedef pointPair<dim> pointPairT;
  typedef double floatT;
  floatT dist;
  pointT u,v;
  pointPair(): dist(floatMax()) {u=pointT();v=pointT();}
  pointPair(pointT uu, pointT vv, floatT distt): u(uu), v(vv), dist(distt) {}
  pointPair(pointT uu, pointT vv): u(uu), v(vv) {
    dist = u.pointDist(v);}
  bool isEmpty() {return dist==floatMax();}
  void closer(pointT uu, pointT vv) {
    floatT distt = uu.pointDist(vv);
    if (distt<dist) {
      u = uu;
      v = vv;
      dist = distt;}
  }
  void closer(pointT uu, pointT vv, floatT distt) {
    if (distt<dist) {
      u = uu;
      v = vv;
      dist = distt;}
  }
  void closer(pointPairT pair) {
    if (pair.dist < dist) {
      u = pair.u;
      v = pair.v;
      dist = pair.dist;}
  }
  friend bool operator<(pointPairT a, pointPairT b) {
    if (a.dist == b.dist) {
      for (int i=0; i<dim; ++i) {
        if(a.u[i]==b.u[i]) continue;
        else return a.u[i]<b.u[i];}
      return false;//
    }
    else return a.dist < b.dist;}
};

//Deprecate 0827
template<int dim>
struct getDist {
  pointPair<dim> *A;
  getDist(pointPair<dim> *AA): A(AA) {}
  double operator() (intT i){
    return A[i].dist;
  }
};

template<int dim>
struct pointPairCmp {
  typedef pointPair<dim> pointPairT;
  bool operator() (pointPairT i, pointPairT j){
    return i.dist < j.dist;
  }
};

// *************************************************************
//    Point manipulation in a grid
// *************************************************************

template<int dim>
inline bool samePoint(double* i, double* j) {
  for (intT ii=0; ii<dim; ++ii) {
    if (abs(i[ii]-j[ii]) > 0.000001) return false;}
  return true;
}

//a less comparator
template<int dim, class pointT>
inline bool pointGridCmp(pointT p1, pointT p2, pointT pMin, floatT r) {
  for(int i=0; i<dim; ++i) {
    intT xx1 = (intT) floor((p1[i]-pMin[i])/r);
    intT xx2 = (intT) floor((p2[i]-pMin[i])/r);
    if (xx1 != xx2) {
      if (xx1 > xx2) return false;
      else return true;}
  }
  return false;
}

//hash function for float array to a cell
template<int dim>
struct hashFloatToCell {
  typedef double floatT;
  typedef point<dim> pointT;
  static const unsigned int prime = -5;
  static const unsigned int mask = -1;
  static const unsigned int range = (1 << 29);
  static const bool noRandom = false;
  int rands[10] = {846930886, 1681692777, 1714636915, 1957747793, 424238335, 719885386, 1649760492, 596516649, 1189641421, 120120309};
  int randInt[dim];
  floatT r;
  pointT pMin;
  hashFloatToCell(pointT pMinn, floatT rr): pMin(pMinn), r(rr) {
    srand(time(NULL));
    for (intT i = 0; i < dim; i++) {
      if(noRandom) randInt[i] = rands[i] % range + 1;
      else randInt[i] = rand() % range + 1;}
  }
  inline uintT primeHash(intT* x, intT n) {
    unsigned long long temp = 0;
    uintT key = 0;
    for (intT i=0; i<n; i++) {
      temp = (long long) x[i] * (long long) randInt[i];
      temp = (temp & mask) + 5 * (temp >> 32);
      if (temp >= prime) temp -= prime;
      temp += key;
      if (temp >= prime) temp -= prime;
      key = (uintT) temp;
    }
    return key;}
  //+1: later cell, -1 earlier cell, 0 same cell
  inline int compareCell(floatT* x1, floatT* x2) {
    for(int i=0; i<dim; ++i) {
      intT xx1 = (intT) floor((x1[i]-pMin[i])/r);
      intT xx2 = (intT) floor((x2[i]-pMin[i])/r);
      if (xx1 != xx2) {
        if (xx1 > xx2) return 1;
        else return -1;
      }}
    return 0;}
  //+1: larger point, -1 smaller point, 0 same point
  inline int comparePoint(floatT* x1, floatT* x2) {
    for(int i=0; i<dim; ++i) {
      if (x1[i] != x2[i]) {
        if (x1[i] > x2[i]) return 1;
        else return -1;
      }}
    return 0;}
  inline uintT hash(floatT *x) {
    intT xx[dim];
    for(int i=0; i<dim; ++i) {
      xx[i] = (intT) floor((x[i]-pMin[i])/r);}
    return primeHash(xx, dim);}
};

//hashStruct for float array object supporting ->coordinate, ../common/ndHash.h
template<int dim, class objT>
struct aFloatHash {
  typedef double floatT;
  typedef hashFloatToCell<dim> hashFunc;
  typedef objT* eType;
  typedef objT* kType;
  hashFunc* hashF;
  objT* e;
  aFloatHash(hashFunc* hashFF):hashF(hashFF) {
    e = new objT();}
  ~aFloatHash() {}
  eType empty() {return e;}
  kType getKey(eType v) {return v;}
  uintT hash(kType c) {
    return hashF->hash(c->coordinate());
  }
  int cmp(kType c1, kType c2) {
    if (c1->isEmpty() || c2->isEmpty()) return 1;
    return hashF->comparePoint(c1->coordinate(), c2->coordinate());
  }
  //inline int diffPoint(floatT* p1, floatT* p2) {return hashF->comparePoint(p1, p2);}
  bool replaceQ(eType c1, eType c2) {return 1;}
};

// *************************************************************
//    Bruteforce
// *************************************************************

template<int dim>
pointPair<dim> bruteForceParallel(point<dim>* P, intT n) {
  static const intT intMax = numeric_limits<intT>::max();
  pointPair<dim>* A = newA(pointPair<dim>, n*n);
  par_for(intT i = 0; i < n; ++ i) {
    for(intT j = 0; j < n; ++ j) {
      if (i == j) {
        A[i+j*n] = pointPair<dim>(P[i], P[j], intMax);
      } else {
        A[i+j*n] = pointPair<dim>(P[i], P[j]);
      }
    }
  }
  intT I = sequence::minIndex<double, intT, getDist<dim>>(0, n*n, getDist<dim>(A));
  auto R = A[I];
  free(A);
  return R;
}

template<int dim>
pointPair<dim> bruteForceSerial(point<dim>* P, intT n) {
  static const intT intMax = numeric_limits<intT>::max();
  pointPair<dim> A = pointPair<dim>(P[0], P[1], intMax);
  par_for(intT i = 0; i < n; ++ i) {
    for(intT j = i+1; j < n; ++ j) {
      A.closer(P[i], P[j]);
    }
  }
  return A;
}

// *************************************************************
//   Misc
// *************************************************************
template<int dim>
point<dim> pMinSerial(point<dim>* items, intT n) {
  point<dim> pMin = point<dim>(items[0].x);
  for(intT p=0; p<n; p++) {
    pMin.minCoords(items[p].x);}
  return pMin;
}

template<int dim>
point<dim> pMinParallel(point<dim>* items, intT n) {
  point<dim> pMin = point<dim>(items[0].x);
  intT P = getWorkers()*8;
  intT blockSize = (n+P-1)/P;
  point<dim> localMin[P];
  for (intT i=0; i<P; ++i) {
    localMin[i] = point<dim>(items[0].x);}
  par_for(intT p=0; p<P; p++) {
    intT s = p*blockSize;
    intT e = min((p+1)*blockSize,n);
    for (intT j=s; j<e; ++j) {
      localMin[p].minCoords(items[j].x);}
  }
  pMin = point<dim>(items[0].x);
  for(intT p=0; p<P; ++p) {
    pMin.minCoords(localMin[p].x);}
  return pMin;
}
#endif
