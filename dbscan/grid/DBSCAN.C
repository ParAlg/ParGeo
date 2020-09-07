// This code is part of the project "Theoretically Efficient and Practical
// Parallel DBSCAN"
// Copyright (c) 2020 Yiqiu Wang, Yan Gu, Julian Shun
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

#include "DBSCAN.h"

#include "shared.h" //todo
#include "grid.h"
#include "cell.h"
#include "pbbs/gettime.h"

using namespace std;

template<int dim>
struct iPoint {
  typedef point<dim> pointT;
  typedef iPoint<dim> iPointT;
  pointT p;
  intT i;
  iPoint(pointT pp, intT ii):p(pp), i(ii) {}
  inline floatT* coordinate() {return p.coordinate();}
  inline floatT coordinate(intT i) {return p.coordinate(i);}
  floatT operator[](int i) {return coordinate(i);}
  inline floatT pointDist(iPointT p2) {return p.pointDist(p2.p);}
  inline bool isEmpty() {return p.isEmpty();}
  inline void setEmpty() {return p.setEmpty();}
};

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the DBSCAN P using a grid algorithm.
 * @param P a point array.
 * @param n length of P.
 * @param epsilon parameter.
 * @param minPts parameter.
 * @return the cluster id of every point (-1: noise, tie-break border points by closest)
 */
template<int dim>
intT* DBSCAN(point<dim>* P, intT n, floatT epsilon, intT minPts) {
  static const bool serial = false;

  typedef point<dim> pointT;
  typedef iPoint<dim> iPointT;
  typedef cell<dim, iPointT> cellT;
  typedef grid<dim, iPointT> gridT;

  cout << "DBSCAN of " << n << ", dim " << dim << " points" << endl;
  cout << "epsilon = " << epsilon << endl;
  cout << "minPts = " << minPts << endl;
  if (n < 2) abort();

  timing t0; t0.start();

  auto PP = newA(iPointT, n);
  par_for(intT i=0; i<n; ++i) {
    PP[i] = iPointT(P[i], i);}

  pointT pMin;
  if (serial) pMin = pMinSerial(P, n);
  else pMin = pMinParallel(P, n);
  auto G = new gridT(n+1, pMin, epsilon/sqrt(dim));
  if(serial) G->insertSerial(PP, n);
  else G->insertParallel(PP, n);//todo, buffer

  cout << "compute-grid = " << t0.next() << endl;

  //mark core
  intT* coreFlag = newA(intT, n);
  par_for(intT i=0; i<n; ++i) coreFlag[i] = -1;

  auto isCore = [&](iPointT p) {
                  coreFlag[p.i] = 1;
                };
  auto fullBox = [&](cellT* c) {
                   if (c->size() >= minPts) {
                     c->pointMap(isCore);}
                 };
  G->allCellMap(fullBox);

  par_for(intT i=0; i<n; ++i) {
    if (coreFlag[i] < 0) {
      intT count = 0;
      auto isCore = [&] (iPointT p) {
                      if(count >= minPts) return true;
                      if(p.pointDist(PP[i]) <= epsilon) {//todo sqrt opt
                        count ++;}
                      return false;};
      G->nghPointMap(PP[i].coordinate(), isCore);
      if (count >= minPts) coreFlag[i] = 1;
      else coreFlag[i] = 0;
    }
  }

  cout << "mark-core = " << t0.next() << endl;

  //cluster core
  intT* cluster = newA(intT, n);
  par_for(intT i=0; i<n; ++i) cluster[i] = -1;

  free(coreFlag);
  delete G;
  free(PP);
  intT* dummy;
  return dummy;
}

template intT* DBSCAN<2>(point<2>*, intT, floatT, intT);
// template intT* DBSCAN<3>(point<3>*, intT, floatT, intT);
// template intT* DBSCAN<4>(point<4>*, intT, floatT, intT);
template intT* DBSCAN<5>(point<5>*, intT, floatT, intT);
// template intT* DBSCAN<6>(point<6>*, intT, floatT, intT);
// template intT* DBSCAN<7>(point<7>*, intT, floatT, intT);
// template intT* DBSCAN<8>(point<8>*, intT, floatT, intT);
// template intT* DBSCAN<9>(point<9>*, intT, floatT, intT);
