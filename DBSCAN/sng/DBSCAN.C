// This code is the C++ version of "Faster DBSCAN via subsampled similarity
// queries" by Heinrich Jiang, Jennifer Jang, Jakub Lacki
//
// Code copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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
#include "pbbs/gettime.h"
#include "pbbs/unionFind.h"
#include "pbbs/utils.h"
#include "pbbs/parallel.h"

template<int dim, class pointT>
intT* coreSNG_BF(pointT* P, intT n, floatT epsilon, intT minPts, floatT s) {
  intT* coreFlag = newA(intT, n);
  intT ns = ceil(n*s);

  par_for (intT i=0; i<n; ++i) {
    coreFlag[i] = 0;
    for (intT j=0; j<ns; ++j) {
      if (P[i].dist(P[utils::hash(i*n+j)%ns]) <= epsilon) {
        coreFlag[i] ++;
      }}
  }

  intT numCore = 0;
  for (intT i=0; i<n; ++i) {
    if (coreFlag[i] >= minPts) {
      coreFlag[i] = 1;
      numCore ++;
    } else coreFlag[i] = 0;
    //cout << coreFlag[i] << " ";
  }
  //cout << endl << endl;
  cout << "bf-num-core = " << numCore << endl;
  return coreFlag;
}


template<int dim, class pointT>
intT* clusterCoreSNG_BF(pointT* P, intT n, floatT epsilon, intT minPts, floatT s, intT* coreFlag) {
  intT ns = ceil(n*s);
  auto uf = unionFind(n);
  par_for (intT i=0; i<n; ++i) {
    for (intT j=0; j<ns; ++j) {
      intT jj = utils::hash(i*n+j)%ns;
      if (coreFlag[i] && coreFlag[jj] && P[i].dist(P[jj]) <= epsilon) {
        uf.link(i,jj);
      }}
  }
  intT* cluster = newA(intT, n);
  par_for (intT i=0; i<n; ++i) {
    auto pi = P[i];
    cluster[i] = -1;
    if(coreFlag[i]) cluster[i] = uf.find(i);
  }
  // for (intT i=0; i<n; ++i) cout << cluster[i] << " ";
  // cout << endl;
  return cluster;
}

template<int dim, class pointT>
void clusterBorderSNG_BF(pointT* P, intT n, floatT epsilon, intT minPts, floatT s, intT* coreFlag, intT* clusterb) {
  intT ns = ceil(n*s);
  floatT thresh = epsilon*epsilon;
  par_for(intT i=0; i<n; ++i) {
    if (!coreFlag[i]) {
      intT cid = -1;
      floatT cDistSqr = floatMax();
      for(intT j=0; j<ns; ++j) {
        intT jj = utils::hash(i*n+j)%ns;
        if (coreFlag[jj]) {
          auto dist = P[i].distSqr(P[jj]);
          if (dist <= thresh && dist < cDistSqr) {
            cid = clusterb[jj];
            cDistSqr = dist;
          }
        }
      }
      clusterb[i] = cid;
    }
  }
}

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the SNG-DBSCAN P using a grid algorithm.
 * @param P a point array.
 * @param n length of P.
 * @param epsilon parameter.
 * @param minPts parameter.
 * @return the cluster id of every point (-1: noise, tie-break border points by closest)
 */
template<int dim>
intT* DBSCAN(point<dim>* P, intT n, floatT epsilon, intT minPts) {
  typedef point<dim> pointT;

  static const floatT s = 0.001;//parameter

  cout << "SNG-DBSCAN of " << n << ", dim " << dim << " points" << endl;
  cout << "epsilon = " << epsilon << endl;
  cout << "minPts = " << minPts << endl;
  cout << "sampling-rate = " << s << endl;
  if (n < 2) abort();

  timing t0; t0.start();

  floatT epsSqr = epsilon*epsilon;
  pointT pMin = pMinParallel(P, n);

  auto coreFlag = coreSNG_BF<dim, pointT>(P, n, epsilon, minPts, s);
  cout << "mark-core-time = " << t0.next() << endl;

  auto cluster = clusterCoreSNG_BF<dim, pointT>(P, n, epsilon, minPts, s, coreFlag);
  cout << "cluster-core-time = " << t0.next() << endl;

  clusterBorderSNG_BF<dim, pointT>(P, n, epsilon, minPts, s, coreFlag, cluster);
  cout << "cluster-border-time = " << t0.next() << endl;

  free(coreFlag);
  return cluster;
}

template intT* DBSCAN<2>(point<2>*, intT, floatT, intT);
template intT* DBSCAN<3>(point<3>*, intT, floatT, intT);
template intT* DBSCAN<4>(point<4>*, intT, floatT, intT);
template intT* DBSCAN<5>(point<5>*, intT, floatT, intT);
template intT* DBSCAN<6>(point<6>*, intT, floatT, intT);
template intT* DBSCAN<7>(point<7>*, intT, floatT, intT);
template intT* DBSCAN<8>(point<8>*, intT, floatT, intT);
template intT* DBSCAN<9>(point<9>*, intT, floatT, intT);
