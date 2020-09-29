// This code is part of the project "Fast Parallel Algorithms for
// Euclidean Minimum SpanningTree and Hierarchical Spatial Clustering"
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

#include "emst.h"
#include "kdTree.h"
#include "wspd.h"
#include "wspdFilter.h"
#include "mark.h"
#include "parallelKruskal.h"
#include "pbbs/gettime.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the EMST P using the filter kruskal algorithm.
 * @param P a point array.
 * @param n length of P.
 * @return a weighted edge array of size n-1
 */
template<int dim>
wEdge<point<dim>>* emst(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  typedef wsp<nodeT> pairT;
  typedef struct nodeT::bcp bcpT;
  //static const bool serial = false;
  cout << "EMST of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;
  t0.start();
  bool paraTree = true;
  treeT* tree = new treeT(P, n, paraTree, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  floatT rhoLo = -0.1;
  floatT beta = 2;
  intT edgesAdded = 0;
  intT numEdges = 0;

  edgeUnionFind *uf = new edgeUnionFind(n);

  floatT wspdTime = 0;
  floatT kruskalTime = 0;
  floatT markTime = 0;

  t0.stop();

  while (edgesAdded < n-1) {
    t0.start();
    // auto *out = new vector<bcpT>();
    // floatT tmp = filterWspdSerial<treeT, nodeT>(beta, rhoLo, out, tree, uf);
    auto output = filterWspdParallel<treeT, nodeT>(beta, rhoLo, tree, uf);
    floatT tmp = output.first;//rho hi
    auto out = output.second;//bcp vector
    wspdTime += t0.stop();

    cout << "---" << endl;
    cout << "beta = " << beta << endl;
    cout << "rho = " << rhoLo << " -- " << tmp << endl;

    numEdges += out->size();

    if (out->size() <= 0) {
      delete out;
      beta *= 2;
      rhoLo = tmp;
      continue;}

    intT edgeCount = out->size();
    cout << "edges = " << edgeCount << endl;

    struct nodeWrap {
      bcpT *E;
      pointT *s;
      nodeWrap(bcpT *Ein, pointT* ss) : E(Ein), s(ss) {}
      //index from offset
      inline intT GetU(intT i) {return E[i].u-s;}
      inline intT GetV(intT i) {return E[i].v-s;}
      inline double GetWeight(intT i) {return E[i].dist;}
    };

    t0.start();
    edgesAdded += parKruskal::mst(nodeWrap(&out->at(0), P), n, out->size(), uf);
    cout << "mst edges = " << edgesAdded << endl;
    kruskalTime += t0.stop();

    t0.start();
    mark(tree->rootNode(), uf, P);
    markTime += t0.stop();

    delete out;
    beta *= 2;
    rhoLo = tmp;
  }

  cout << endl << "wspd-time = " << wspdTime << endl;
  cout << "kruskal-time = " << kruskalTime << endl;
  cout << "mark-time = " << markTime << endl;

  typedef wEdge<point<dim>> outT;
  auto R = newA(outT, n-1);
  par_for(intT i=0; i<n-1; ++i) {
    auto e = uf->getEdge(i);
    R[i] = outT(P[e.first], P[e.second], P[e.first].dist(P[e.second]));
  }
  cout << "copy-time = " << t0.stop() << endl;

  delete uf;
  return R;
}

template wEdge<point<2>>* emst(point<2>*, intT);
template wEdge<point<3>>* emst(point<3>*, intT);
template wEdge<point<4>>* emst(point<4>*, intT);
template wEdge<point<5>>* emst(point<5>*, intT);
template wEdge<point<6>>* emst(point<6>*, intT);
template wEdge<point<7>>* emst(point<7>*, intT);
template wEdge<point<8>>* emst(point<8>*, intT);
template wEdge<point<9>>* emst(point<9>*, intT);
