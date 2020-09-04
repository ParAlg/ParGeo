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

#include "sieveGrid.h"
#include "pbbs/randPerm.h"
using namespace std;

// *************************************************************
//    Rabin's algorithm
// *************************************************************

/**
 * Computes the closest pair of P using Rabin's randomized algorithm.
 * @param PIn a point array.
 * @param n length of PIn.
 * @param permuted set to true if PIn is already permuted.
 * @param serial optional boolean variable which should be set to true if wish to run serialy.
 * @return the closest pair
 */
template<int dim>
pointPair<dim> rabin(point<dim>* PIn, intT n, bool permuted=false, bool serial=false) {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointPair<dim> pointPairT;
  static const bool noRandom = false;
  static const bool verbose = true;
  static const double xx0 = 0.8;
  bool useTree = false;
  if(dim>=5) useTree = true;//turns on kdTree when dim>=5
  if(n<1000) return bruteForceSerial<dim>(PIn, n);

  pointT pMin;
  if (serial) {
    pMin = pMinSerial(PIn, n);
  } else {
    pMin = pMinParallel(PIn, n);
  }

  if(!noRandom && !permuted) {
    //permutation
    if(verbose) cout << "permuting points" << endl;
    randPerm(PIn, n);
  }

  intT x0 = intT(pow(n, xx0));
  cout << "sample size = " << x0 << endl;
  auto R0 = rabin(PIn, x0, true, serial);//recursively call rabin's algorithm

  auto grids = grid<dim>(n, pMin, R0.dist);
  auto hash = cellHash<dim>(&grids);
  Table<cellHash<dim>,intT>* table;
  if (useTree) table = new Table<cellHash<dim>,intT>(1, hash);
  else table = new Table<cellHash<dim>,intT>(n, hash, 4);
  grids.addTable(table);
  if(verbose) cout << "r = " << R0.dist << endl;

  if(serial) grids.reInsertSerial(PIn, n, useTree);
  else grids.reInsertParallel(PIn, n, NULL, useTree);
  if(serial) return grids.closestPairSerial(useTree);
  else return grids.closestPairParallel(useTree);
}

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the closest pair of P using Rabin's randomized algorithm.
 * @param P a point array.
 * @param n length of P.
 * @return the closest pair
 */
template<int dim>
pair<point<dim>, point<dim>> closestPair(point<dim>* P, intT n) {
  static const bool serial = false;
  cout << "closestPair of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  pointPair<dim> r = rabin<dim>(P, n, false, serial);
  cout << "rabin's " << r.u << ", " << r.v << ", dist " << r.dist << endl << endl;

  return make_pair(r.u, r.v);
}

template pair<point<2>, point<2>> closestPair<2>(point<2>*, intT);
template pair<point<3>, point<3>> closestPair<3>(point<3>*, intT);
template pair<point<4>, point<4>> closestPair<4>(point<4>*, intT);
template pair<point<5>, point<5>> closestPair<5>(point<5>*, intT);
template pair<point<6>, point<6>> closestPair<6>(point<6>*, intT);
template pair<point<7>, point<7>> closestPair<7>(point<7>*, intT);
template pair<point<8>, point<8>> closestPair<8>(point<8>*, intT);
template pair<point<9>, point<9>> closestPair<9>(point<9>*, intT);
