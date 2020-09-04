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
#include "pbbs/gettime.h"
using namespace std;

// *************************************************************
//    Sieve algorithm
// *************************************************************
static constexpr bool noRandom = false;
intT rands[10] = {1,111,1100,1012,300123,1230123,12031230,1231239,5430,12093812};

//computes the closest neighbor distance of a random point in P
template<int dim, class floatT>
floatT getDiSerial(point<dim>* P, intT nn) {
  intT xi;
  if(noRandom) xi = rands[2]%nn;
  else xi = rand()%nn;
  floatT r = floatMax();
  for(intT i=0; i<nn; ++i) {
    if (i == xi) continue;
    floatT dist = P[i].pointDist(P[xi]);
    r = r > dist ? dist : r;}
  return r;}

//computes the closest neighbor distance of a random point in P
template<int dim, class floatT>
floatT getDiParallel(point<dim>* P, intT nn, floatT* xDists) {
  intT xi;
  if(noRandom) xi = rands[2]%nn;
  else xi = rand()%nn;
  par_for(intT i=0; i<nn; ++i) {
    xDists[i] = P[i].pointDist(P[xi]);
  }
  xDists[xi] = xDists[0]>xDists[1] ? xDists[0] : xDists[1];//make sure xi's dist to self (zero) isn't counted
  auto getItem = [&](intT i) {return xDists[i];};
  intT iMin = sequence::minIndex<floatT,intT>(0, nn, getItem);
  return xDists[iMin];}

/**
 * Computes the closest pair of P using the sieve randomized algorithm.
 * @param PIn a point array.
 * @param n length of PIn.
 * @param serial whether to run serially.
 * @return the closest pair
 */
template<int dim>
pointPair<dim> sieve(point<dim>* PIn, intT n, bool serial=false) {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointPair<dim> pointPairT;

  static const bool verbose = true;
  bool useTree = false;
  if(dim>=5) useTree = true;

  //make a copy of PIn in P, where P will be modified
  pointT* P = newA(pointT, n);
  par_for(intT i=0; i<n; ++i) {
    P[i] = PIn[i];
  }

  timing t0; t0.start();
  double tGetMin = 0;
  double tInsert = 0;
  double tCheck = 0;

  floatT* xDists = newA(floatT, n);
  pointT* PP = newA(pointT, n);
  intT nn = n;

  if(!noRandom) {
    //permutation
    if(verbose) cout << "permuting points" << endl;
    randPerm(P, n);
  }

  //compute pMin and initial grid size
  pointT pMin;
  floatT r;
  if (serial) {
    pMin = pMinSerial(P, nn);
    r = getDiSerial<dim,floatT>(P, nn);
  } else {
    pMin = pMinParallel(P, nn);
    r = getDiParallel(P, nn, xDists);
  }
  r /= 3;

  auto grids = grid<dim>(n, pMin, r);
  auto hash = cellHash<dim>(&grids);
  Table<cellHash<dim>,intT>* table;
  if (useTree) table = new Table<cellHash<dim>,intT>(1, hash);
  else table = new Table<cellHash<dim>,intT>(n, hash, 4);
  grids.addTable(table);
  if(verbose) cout << "r = " << r << endl;
  cout << "init-time = " << t0.stop() << endl;

  intT* aux0 = newA(intT, n+1);
  bool* aux1 = newA(bool, n+1);
  intT rounds = 0;
  while (nn > 2) {
    rounds++;
    t0.start();
    if (serial) r = getDiSerial<dim,floatT>(P, nn);
    else r = getDiParallel(P, nn, xDists);
    r /= 3;
    if(verbose) cout << "r = " << r << ", ";
    grids.updateR(r);
    tGetMin += t0.next();
    if (serial) grids.reInsertSerial(P, nn, useTree);
    else grids.reInsertParallel(P, nn, aux0, useTree);
    tInsert += t0.next();
    _seq<pointT> pSeq;
    if (serial) pSeq = grids.checkInsertSerial(P, PP, nn, useTree);
    else pSeq = grids.checkInsertParallel(P, PP, nn, aux1, useTree);
    nn = pSeq.n;
    if(verbose) cout << "rem-points = " << nn << endl;
    swap(P,PP);
    tCheck += t0.stop();
  }

  if(verbose) cout << "-- final insert" << endl;
  r*=3;
  if(verbose) cout << "r = " << r << endl;
  grids.updateR(r);
  t0.start();
  if (serial) grids.reInsertSerial(PIn, n, useTree);
  else grids.reInsertParallel(PIn, n, aux0, useTree);
  tInsert += t0.stop();
  pointPairT R;
  if (serial) R = grids.closestPairSerial(useTree);
  else R = grids.closestPairParallel(useTree);

  cout << "getmin-time = " << tGetMin << endl;
  cout << "insert-time = " << tInsert << endl;
  cout << "check-time = " << tCheck << endl;
  cout << rounds << " rounds" << endl;

  free(aux0);
  free(aux1);
  free(xDists);
  free(P);
  return R;
}

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the closest pair of P using the sieve randomized algorithm.
 * @param P a point array.
 * @param n length of P.
 * @return the closest pair
 */
template<int dim>
pair<point<dim>, point<dim>> closestPair(point<dim>* P, intT n) {
  static const bool serial = false;
  cout << "closestPair of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  pointPair<dim> r = sieve<dim>(P, n, serial);
  cout << "sieve " << r.u << ", " << r.v << ", dist " << r.dist << endl << endl;

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
