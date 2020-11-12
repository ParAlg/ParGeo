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
#include "pbbs/unionFind.h"
#include "pbbs/utils.h"
#include "pbbs/parallel.h"
#include "kdTree.h"
#include "kdNode.h"
#include "coreBccp.h"
#include "bruteforce.h"

//todo in-place kdtree (if slow)
//todo cluster core slow for large eps todo check again
//todo jemalloc

using namespace std;

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
  static const bool checker = false;

  typedef point<dim> pointT;
  typedef cell<dim, pointT> cellT;
  typedef grid<dim, pointT> gridT;

  cout << "DBSCAN of " << n << ", dim " << dim << " points" << endl;
  cout << "epsilon = " << epsilon << endl;
  cout << "minPts = " << minPts << endl;
  if (n < 2) abort();

  timing t0; t0.start();

  floatT epsSqr = epsilon*epsilon;
  pointT pMin = pMinParallel(P, n);

  auto G = new gridT(n+1, pMin, epsilon/sqrt(dim));
  G->insertParallel(P, n);
  cout << "num-cell = " << G->numCell() << endl;
  cout << "compute-grid = " << t0.next() << endl;

  //mark core
  intT* coreFlag = newA(intT, n);
  par_for(intT i=0; i<n; ++i) coreFlag[i] = -1;

  auto isCore = [&](pointT *p) {
                  coreFlag[p-P] = 1;
                  return false;
                };

  par_for(intT i=0; i<G->numCell(); ++i) {
    cellT* c = G->getCell(i);
    if (c->size() >= minPts) c->pointMap(isCore);
  }

  par_for(intT i=0; i<n; ++i) {
    if (coreFlag[i] < 0) {
      intT count = 0;
      auto isCore = [&] (pointT *p) {
                      if(count >= minPts) return true;
                      if(p->distSqr(P[i]) <= epsSqr) {//todo sqrt opt
                        count ++;}
                      return false;};
      G->nghPointMap(P[i].coordinate(), isCore);
      if (count >= minPts) coreFlag[i] = 1;
      else coreFlag[i] = 0;
    }
  }

  cout << "mark-core-time = " << t0.next() << endl;

  if (checker) {
    auto cfCheck = coreBF<dim, pointT>(P, n, epsilon, minPts);
    intT numCore = 0;
    for (intT i=0; i<n; ++i) {
      if (cfCheck[i] != coreFlag[i]) {
        cout << "error, core flag mismatch, abort()" << endl;
        abort();
      }
      if(coreFlag[i]) numCore ++;
    }
    cout << "num-core = " << numCore << endl;
    free(cfCheck);
    cout << "mark-core-correctness-checked = " << t0.next() << endl;
  }

  //cluster core
  auto ccFlag = newA(intT, G->numCell());
  par_for(intT i=0; i<G->numCell(); ++i) {
    auto ci = G->getCell(i);
    ccFlag[i] = 0;
    auto hasCore = [&](pointT *p) {
                     if (coreFlag[p-P]) {
                       ccFlag[i] = 1;
                       return true;
                     }
                     return false;
                   };
    ci->pointMap(hasCore);
  }

  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  typedef typename nodeT::bcp bcpT;
  auto trees = newA(treeT*, G->numCell());
  par_for(intT i=0; i<G->numCell(); ++i) trees[i] = NULL;

  // auto degCmp = [&](intT i, intT j) {
  //                 return G->getCell(i)->size() < G->getCell(j)->size();
  //               };
  // auto ordering = newA(intT, G->numCell());
  // par_for(intT i=0; i<G->numCell(); ++i) ordering[i] = i;
  //sampleSort(ordering, G->numCell(), degCmp);

  auto uf = unionFind(G->numCell());

  floatT bcpTotalTime = 0; timing t1;
  par_for(intT i=0; i<G->numCell(); ++i) {
    if (ccFlag[i]) {
      auto ti = trees[i];
      auto procTj = [&](cellT* cj) {
                      intT j = cj - G->getCell(0);
                      auto tj = trees[j];
                      if (j < i && ccFlag[j] &&
                          uf.find(i) != uf.find(j)) {
                        if(hasEdge<cellT, treeT, pointT>(i, j, coreFlag, P, epsilon, G->getCell(0), trees)) {
                          uf.link(i, j);
                        }
                      }
                      return false;
                    };
      //G->nghCellMap(G->getCell(ordering[i]), procTj);
      G->nghCellMap(G->getCell(i), procTj);
    }
  }

  par_for(intT i=0; i<G->numCell(); ++i) {
    if (trees[i]) delete trees[i];
  }

  intT* cluster = newA(intT, n);
  par_for(intT i=0; i<n; ++i) cluster[i] = -1;

  par_for(intT i=0; i<G->numCell(); ++i) {
    auto cid = G->getCell(uf.find(i))->getItem() - P;//id of first point
    auto clusterCore = [&](pointT* p){
                         if (coreFlag[p - P])
                           cluster[p - P] = cid;
                         return false;
                       };
    G->getCell(i)->pointMap(clusterCore);
  }
  cout << "cluster-core-time = " << t0.next() << endl;

  if (checker) {
    auto ccCheck = clusterCoreBF<dim, pointT>(P, n, epsilon, minPts, coreFlag);
    for (intT i=0; i<n; ++i) {
      for (intT j=i+1; j<n; ++j) {
        bool sameClu = (cluster[i] >= 0) && (cluster[i] == cluster[j]);
        bool sameCluCheck = (ccCheck[i] >= 0) && (ccCheck[i] == ccCheck[j]);
        if (sameClu != sameCluCheck) {
          cout << "error, cluster core mismatch, abort()" << endl;
          abort();
        }
      }
    }
    free(ccCheck);
    cout << "cluster-core-correctness-checked = " << t0.next() << endl;
  }

  intT* cbCheck; if (checker) {
    cbCheck = newA(intT, n);
    par_for(intT i=0; i<n; ++i) cbCheck[i] = cluster[i];
  }

  //cluster border to closest core point
  par_for(intT i=0; i<n; ++i) {
    if (!coreFlag[i]) {
      intT cid = -1;
      floatT cDistSqr = floatMax();
      auto closestCore = [&] (pointT* p) {
                           if (coreFlag[p-P]) {
                             auto dist = p->distSqr(P[i]);
                             if (dist <= epsSqr && dist < cDistSqr) {
                               cDistSqr = dist;
                               cid = cluster[p-P];}
                           }
                           return false;};
      G->nghPointMap(P[i].coordinate(), closestCore);
      cluster[i] = cid;
    }
  }
  //for (intT i=0; i<n; ++i) cout << cluster[i] << " ";cout << endl << endl;
  cout << "cluster-border-time = " << t0.next() << endl;

  if (checker) {
    clusterBorderBF<dim, pointT>(P, n, epsilon, minPts, coreFlag, cbCheck);
    //for (intT i=0; i<n; ++i) cout << cbCheck[i] << " ";cout << endl << endl;
    for (intT i=0; i<n; ++i) {
      for (intT j=i+1; j<n; ++j) {
        bool sameClu = (cluster[i] >= 0) && (cluster[i] == cluster[j]);
        bool sameCluCheck = (cbCheck[i] >= 0) && (cbCheck[i] == cbCheck[j]);
        if (sameClu != sameCluCheck) {
          cout << "error, cluster border mismatch, abort()" << endl;
          abort();
        }
      }
    }
    free(cbCheck);
    cout << "cluster-border-correctness-checked = " << t0.next() << endl;
  }

  free(coreFlag);
  free(ccFlag);
  free(trees);
  delete G;
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
