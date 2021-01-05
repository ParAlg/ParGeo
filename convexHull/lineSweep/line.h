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
    I = newA(intT, n+1);//note +1
  }

  timing t; t.start();

  auto yCmp = [&](point2d& a, point2d& b) {
    return a.y() < b.y();
  };
  sampleSort(&P[0], n, yCmp);
  cout << "sort-time = " << t.next() << endl;

  floatT xLeft = min(P[0].x(), P[n-1].x());
  floatT xRight = max(P[0].x(), P[n-1].x());
  cout << "x-left = " << xLeft << endl;
  cout << "x-right = " << xRight << endl;

  intT ml = 0;
  intT mr = n;

  auto pushLeft = [&](intT i) {I[ml++] = i;};
  auto popLeft = [&]() {ml--;};
  auto pushRight = [&](intT i) {I[mr--] = i;};
  auto popRight = [&]() {mr++;};
  auto sizeLeft = [&]{return ml;};
  auto sizeRight = [&]{return n-mr;};
  pushLeft(0);
  pushRight(0);

  //scan from bottom up
  for (intT i=1; i<n; ++i) {
    if (P[i].x()>xLeft && P[i].x()<xRight) {
      continue;
    } else if (P[i].x()<=xLeft) {
      if (sizeLeft() >= 2 && !rightTurn(P[I[ml-2]], P[I[ml-1]], P[i])) {
	popLeft();
	while (sizeLeft()>1 && !rightTurn(P[I[ml-2]], P[I[ml-1]], P[i])) popLeft();
      }
      pushLeft(i);
    } else {//>=xRight
      if (sizeRight() >= 2 && !leftTurn(P[I[mr+2]], P[I[mr+1]], P[i])) {
	popRight();
	while (sizeRight()>1 && !leftTurn(P[I[mr+2]], P[I[mr+1]], P[i])) popRight();
      }
      pushRight(i);
    }
  }

  cout << sizeLeft() + sizeRight() << endl;
  cout << "scan-time = " << t.stop() << endl;

  //todo reorder

  return _seq<intT>(I, 1);//todo
}

#endif
