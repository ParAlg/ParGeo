// Copyright (c) 2020 Yiqiu Wang
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

#include "spanner.h"
#include "kdTree.h"
#include "wspd.h"
#include "wspdNormal.h"
#include "pbbs/gettime.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the Spanner of P using the WSPD.
 * @param P a point array.
 * @param n length of P.
 * @param t spanner parameter t, maximum relaxation factor of path wrt spatial distance.
 * @return a weighted edge array of size n-1
 */
template<int dim>
wEdge<point<dim>>* spanner(point<dim>* P, intT n, intT t) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  typedef wsp<nodeT> pairT;
  typedef struct nodeT::bcp bcpT;

  static const bool serial = true;

  cout << t << "-spanner of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  floatT s = 4*((floatT)t+1)/((floatT)t-1);
  cout << "separation-constant = " << s << endl;
  timing t0;
  t0.start();
  bool paraTree = true;
  treeT* tree = new treeT(P, n, paraTree, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  vector<pairT> *wspd;
  if(serial) {
    wspd = new vector<pairT>();
    auto wg = wspdNormalSerial<nodeT>(wspd);
    wspdSerial<nodeT, wspdNormalSerial<nodeT>>(tree->rootNode(), &wg, s);
  } else {
    auto wgpar = wspdNormalParallel<nodeT>(tree->rootNode()->size());
    wspdParallel<nodeT, wspdNormalParallel<nodeT>>(tree->rootNode(), &wgpar, s);
    wspd = wgpar.collect();
  }

  cout << "#wsp = " << wspd->size() << endl;
  cout << "wspd-time = " << t0.next() << endl;

  typedef wEdge<point<dim>> outT;
  auto R = newA(outT, wspd->size());
  par_for(intT i=0; i<wspd->size(); ++i) {
    auto pt1 = wspd->at(i).u->getItem(0);
    auto pt2 = wspd->at(i).v->getItem(0);
    R[i] = outT(*pt1, *pt2, pt1->dist(*pt2));
  }
  cout << "copy-time = " << t0.stop() << endl;
  cout << "edge-count = " << wspd->size() << endl;

  delete wspd;
  return R;
}

template wEdge<point<2>>* spanner(point<2>*, intT, intT);
template wEdge<point<3>>* spanner(point<3>*, intT, intT);
template wEdge<point<4>>* spanner(point<4>*, intT, intT);
template wEdge<point<5>>* spanner(point<5>*, intT, intT);
template wEdge<point<6>>* spanner(point<6>*, intT, intT);
template wEdge<point<7>>* spanner(point<7>*, intT, intT);
template wEdge<point<8>>* spanner(point<8>*, intT, intT);
template wEdge<point<9>>* spanner(point<9>*, intT, intT);
