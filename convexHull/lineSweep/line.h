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

#ifndef LINE_H
#define LINE_H

#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "hull.h"

_seq<intT> lineSweepSerial(point2d* P, intT n, intT* I=NULL) {
  static const bool verbose = false;

  //whether a,b,c forms a right turn
  auto rightTurn = [&](point2d& a, point2d& b, point2d& c) {
    auto cross = (b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x());
    return cross <= 0;
  };
  auto leftTurn = [&](point2d& a, point2d& b, point2d& c) {
    return !rightTurn(a,b,c);
  };

  if (!I) {
    I = newA(intT, n);
  }

  timing t1;t1.start();

  auto yCmp = [&](point2d& a, point2d& b) {
    return a.y() < b.y();
  };
  sampleSort(&P[0], n, yCmp);
  cout << " sort-time = " << t1.next() << endl;

  intT ml = 0;
  intT mr = n-1;

  auto pushLeft = [&](intT i) {I[ml++] = i;};
  auto popLeft = [&]() {ml--;};
  auto pushRight = [&](intT i) {I[mr--] = i;};
  auto popRight = [&]() {mr++;};
  auto sizeLeft = [&]{return ml;};
  auto sizeRight = [&]{return n-1-mr;};
  pushLeft(0);

  //scan from bottom up
  for (intT i=1; i<n; ++i) {
    if (triArea(P[0], P[n-1], P[i]) > numericKnob) {
      while (sizeLeft() >= 2 && !rightTurn(P[I[ml-2]], P[I[ml-1]], P[i]))
	popLeft();
      pushLeft(i);
    } else {
      while (1) {
	if (sizeRight() >= 2 && !leftTurn(P[I[mr+2]], P[I[mr+1]], P[i]))
	  popRight();
	else if (sizeRight() == 1 && !leftTurn(P[I[0]], P[I[mr+1]], P[i]))
	  popRight();
	else break;
      }
      pushRight(i);
    }
  }

  while (sizeLeft() >= 2 && !rightTurn(P[I[ml-2]], P[I[ml-1]], P[I[mr+1]]))
    popLeft();

  cout << " scan-time = " << t1.next() << endl;

  auto AR = I+n-sizeRight();
  for(intT i=0; i<sizeRight(); ++i) {
    I[sizeLeft()+i] = AR[i];
  }

  cout << " move-time = " << t1.stop() << endl;

  return _seq<intT>(I, sizeLeft()+sizeRight());
}

#endif
