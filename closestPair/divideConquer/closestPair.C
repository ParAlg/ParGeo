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

#include "divideConquer.h"
using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the closest pair of P using a divide and conquer algorithm.
 * @param P a point array.
 * @param n length of P.
 * @return the closest pair
 */
template<int dim>
pair<point<dim>, point<dim>> closestPair(point<dim>* P, intT n) {
  static const bool serial = false;
  cout << "closestPair of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();
  // pointPair<dim> rr = bruteForceSerial<dim>(P, n);
  // cout << "bruteforce " << rr.u << ", " << rr.v << ", dist " << rr.dist << endl << endl;

  point<dim>* Pp = newA(point<dim>, n);
  pointPair<dim> r;
  if (serial) r = divideConquerSerial<dim>(P, n, Pp, dim-1);
  else r = divideConquerParallel<dim>(P, n, Pp, dim-1);
  free(Pp);
  cout << "divideConquer " << r.u << ", " << r.v << ", dist " << r.dist << endl << endl;

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
