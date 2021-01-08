// This code is part of the Pargeo Library
// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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

#include <algorithm>
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "quick.h"
#include "gift.h"
#include "graham.h"
using namespace std;
using namespace sequence;

_seq<intT> hullWorker(point2d* P, intT n, intT* I, intT* newI, floatT& bfTime) {
  typedef _seq<intT> hullT;

  if (n <= 1000) {//todo tune
    auto H = quickHullSerial(P, n, I);
    // cout << "sub n = " << H.n << endl;
    // for(intT i=0; i<H.n; ++i) {
    //   cout << P[I[i]] << " ";
    // }
    // cout << endl;
    return H;
  }

  intT n1 = intT(n/2);

  hullT H1, H2;
  par_do([&](){H1 = hullWorker(P, n1, I, newI+n1, bfTime);},
	 [&](){H2 = hullWorker(P+n1, n-n1, I+n1, newI+n1, bfTime);});
  for(intT i=0; i<H2.n; ++i) H2.A[i] += n1;

  //check if p is above line(l,r)
  auto aboveLine = [&](intT p, intT l, intT r) {
		     return triArea(P[p], P[l], P[r]) > numericKnob;};

  static const int UpTangent = 1;
  static const int LowTangent = 2;
  static const int Left = 3;
  static const int Right = 4;

  //check point a \in H1, wrt (a,b) where b \in H2
  // or point a \in H2, wrt (a,b) where b \in H1
  auto check = [&](intT c1, intT c2, intT a, intT b) {
		 auto b1 = aboveLine(c1, a, b);
		 auto b2 = aboveLine(c2, a, b);
		 // cout << "line = (" << a << ": " << P[a] << ", "  << b << ": " << P[b] << "), ";
		 // cout << "p1 = (" << c1 << ")" << P[c1] << ", ";
		 // cout << "p2 = (" << c2 << ")" << P[c2] << ", ";

		 if (!b1 && !b2) {
		   //cout << "up" << endl;
		   return UpTangent;
		 } else if (!b1 && b2) {
		   //cout << "left" << endl;
		   return Left;
		 } else if (b1 && !b2) {
		   //cout << "right" << endl;
		   return Right;
		 } else {
		   //cout << "low" << endl;
		   return LowTangent;
		 }
	       };

  intT u1, u2, l1, l2;
  // cout << "n = " << n << endl;
  // cout << "size1 = " << H1.n << endl;
  // for(intT i=0; i<H1.n; ++i) cout << H1.A[i] << " "; cout << endl;
  // cout << "size2 = " << H2.n << endl;
  // for(intT i=0; i<H2.n; ++i) cout << H2.A[i] << " "; cout << endl;

  timing t; t.start();
  for(intT i=0; i<H1.n; ++i) {
    for(intT j=0; j<H2.n; ++j) {
      int checkPrev = check(H1.A[(i-1+H1.n)%H1.n],H1.A[(i+1)%H1.n],H1.A[i],H2.A[j]);
      int checkNext = check(H2.A[(j-1+H2.n)%H2.n],H2.A[(j+1)%H2.n],H1.A[i],H2.A[j]);
      if (checkPrev==UpTangent && checkNext==UpTangent) {
	u1 = i;
	u2 = j;
	//cout << "upper found " << u1 << ", " << u2 << endl;
      }
      if (checkPrev==LowTangent && checkNext==LowTangent) {
	l1 = i;
	l2 = j;
	//cout << "lower found " << l1 << ", " << l2 << endl;
      }
    }
  }
  // cout << "up-tangent = " << H1.A[u1] << ", " << H2.A[u2] << endl;
  // cout << "low-tangent = " << H1.A[l1] << ", " << H2.A[l2] << endl;
  bfTime += t.stop();

  intT ii = 0;
  intT i = l1;
  do {
    newI[ii++] = H1.A[i%H1.n];
  } while (i++%H1.n != u1);
  if(newI[ii-1]!=H1.A[u1]) newI[ii++] = H1.A[u1];

  i = u2;
  do {
    newI[ii++] = H2.A[i%H2.n];
  } while (i++%H2.n != l2);
  if(newI[ii-1]!=H2.A[l2]) newI[ii++] = H2.A[l2];

  for(intT i=0; i<ii; ++i) {
    H1.A[i] = newI[i];
    //cout << newI[i] << " ";
  }
  H1.n = ii;
  //cout << endl;

  return H1;
}

_seq<intT> hull(point2d* P, intT n) {

  timing t; t.start();

  auto xLess = [&](point2d p1, point2d p2) {
		 return p1.x()<p2.x();};
  timing t0; t0.start();
  sampleSort(P, n, xLess);
  cout << "sort-time = " << t0.stop() << endl;

  floatT bfTime = 0;
  auto I = newA(intT, n);
  auto I2 = newA(intT, n);
  auto H = hullWorker(P, n, I, I2, bfTime);
  free(I2);

  floatT totalTime = t.stop();
  cout << "bf-merge-time = " << bfTime << "(" << 100*bfTime/totalTime << "%)" << endl;
  cout << "total-hull-time = " << totalTime << endl;

  check(P, n, H.A, H.n);

  return H;
}
