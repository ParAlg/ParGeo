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

#include "hdbscan.h"
#include "kdTree.h"
#include "wspd.h"
#include "unreachableDecomp.h"
#include "mark.h"
#include "parallelKruskal.h"
#include "kNearestNeighbors.h"
#include "kBuffer.h"
#include "dendrogram.h"
#include "pbbs/gettime.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the HDBSCAN* P using the filter kruskal algorithm.
 * @param P a point array.
 * @param n length of P.
 * @return a weighted edge array of size n-1
 */
template<int dim>
dEdge* hdbscan(point<dim>* P, intT n, floatT eps, intT minPts) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  typedef wsp<nodeT> pairT;
  typedef struct nodeT::bcp bcpT;
  //static const bool serial = false;
  cout << "HDBSCAN of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;
  t0.start();
  bool paraTree = true;
  treeT* tree = new treeT(P, n, paraTree, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  floatT* coreDist = newA(floatT, n);
  auto buffer = KBuffer::allocKBuffer<pointT*>(minPts, n);
  parallel_for (0, n, [&](intT i) {
			auto buf = buffer + i;
			tree->rootNode()->kNN(&P[i], minPts, buf);
			coreDist[i] = buf->get(minPts-1)->dist(P[i]);
		      });
  cout << "knn-time = " << t0.next() << endl;

  floatT* cdMin = newA(floatT, tree->size()*2);
  floatT* cdMax = newA(floatT, tree->size()*2);
  parallel_for(0, tree->size()*2,
	       [&](intT i) {
		 cdMin[i] = floatMax();
		 cdMax[i] = floatMin();
	       });
  nodeCD(tree->rootNode(), coreDist, cdMin, cdMax, tree->rootNode(), P);

  floatT rhoLo = -0.1;
  floatT beta = 2;
  intT edgesAdded = 0;
  intT numEdges = 0;

  edgeUnionFind *uf = new edgeUnionFind(n);

  while (edgesAdded < n-1) {
    //auto *out = new vector<bcpT>();
    //floatT tmp = filterUnreachSerial<dim, treeT, nodeT, pointT>(beta, rhoLo, out, tree, uf, coreDist, cdMin, cdMax, P);

    auto output = filterUnreachParallel<dim, treeT, nodeT, pointT>(beta, rhoLo, tree, uf, coreDist, cdMin, cdMax, P);
    floatT tmp = output.first;//rho hi
    auto out = output.second;//bcp vector
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

    edgesAdded += parKruskal::mst(nodeWrap(&out->at(0), P), n, out->size(), uf);
    cout << "mst edges = " << edgesAdded << endl;

    mark(tree->rootNode(), uf, P);

    delete out;
    beta *= 2;
    rhoLo = tmp;
  }

  typedef dEdge outT;
  auto R = newA(outT, n-1);
  parallel_for(0, n-1,
	       [&](intT i) {
		 auto e = uf->getEdge(i);
		 auto w = max(P[e.first].dist(P[e.second]), coreDist[e.first]);
		 w = max(w, coreDist[e.second]);
		 R[i] = outT(e.first, e.second, w);
	       });
  cout << "copy-time = " << t0.stop() << endl;

  auto myDendro = dendrogram::directedDendro<dEdge>(R, n);
  cout << "dendrogram-time = " << t0.stop() << endl;

  free(coreDist);
  free(cdMin);
  free(cdMax);
  delete uf;
  return R;
}

template dEdge* hdbscan(point<2>*, intT, floatT, intT);
template dEdge* hdbscan(point<3>*, intT, floatT, intT);
template dEdge* hdbscan(point<4>*, intT, floatT, intT);
template dEdge* hdbscan(point<5>*, intT, floatT, intT);
template dEdge* hdbscan(point<6>*, intT, floatT, intT);
template dEdge* hdbscan(point<7>*, intT, floatT, intT);
template dEdge* hdbscan(point<8>*, intT, floatT, intT);
template dEdge* hdbscan(point<9>*, intT, floatT, intT);
