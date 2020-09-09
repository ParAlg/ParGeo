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

#include "emst.h"
#include "kdTree.h"
#include "wspd.h"
#include "wspdNormal.h"
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

  auto wgpar = wspdNormalParallel<nodeT>(tree->rootNode()->size());
  wspdParallel<nodeT, wspdNormalParallel<nodeT>>(tree->rootNode(), &wgpar);
  auto out = wgpar.collect();
  // vector<pairT> *out = new vector<pairT>();
  // auto wg = wspdNormalSerial<nodeT>(out);
  // wspdSerial<nodeT, wspdNormalSerial<nodeT>>(tree->rootNode(), &wg);

  cout << "#wsp = " << out->size() << endl;
  cout << "wspd-time = " << t0.next() << endl;

  struct indexBcp {
    intT u,v;
    floatT dist;
    indexBcp(intT uu, intT vv, floatT distt): u(uu), v(vv), dist(distt) {};
  };

  auto bcps = newA(indexBcp, out->size());

  par_for (intT i=0; i<out->size(); ++i) {
    bcpT bcp = out->at(i).u->compBcp(out->at(i).v);
    bcps[i] = indexBcp(bcp.u-P, bcp.v-P, bcp.dist);
  }
  cout << "bcp-time = " << t0.next() << endl;

  edgeUnionFind *uf = new edgeUnionFind(n);
  parKruskal::mst(parKruskal::bcpWrap<indexBcp>(bcps), n, out->size(), uf);
  cout << "kruskal-time = " << t0.next() << endl;

  typedef wEdge<point<dim>> outT;
  auto R = newA(outT, n-1);
  par_for(intT i=0; i<n-1; ++i) {
    auto e = uf->getEdge(i);
    R[i] = outT(P[e.first], P[e.second], P[e.first].dist(P[e.second]));
  }
  cout << "copy-time = " << t0.stop() << endl;

  delete uf;
  free(bcps);
  delete out;
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
