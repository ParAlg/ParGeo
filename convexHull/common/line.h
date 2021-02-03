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

namespace lineInternal {

  //whether a,b,c forms a right/left turn, disallows straight line
  inline bool rightTurn(point2d& a, point2d& b, point2d& c) {
    auto cross = (b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x());
    return cross < 0;
  };

  inline bool leftTurn(point2d& a, point2d& b, point2d& c) {
    auto cross = (b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x());
    return cross > 0;
  };

}

intT lineSweepSerial(point2d* P, intT n, intT* I) {
  static const bool verbose = false;

  if (!I) abort();

#ifndef SILENT
  timing t1;t1.start();
#endif

  auto yCmp = [&](point2d& a, point2d& b) {
    return a.y() < b.y();
  };
  quickSortSerial(&P[0], n, yCmp);
#ifndef SILENT
  cout << " sort-time = " << t1.next() << endl;
#endif

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
      while (sizeLeft() >= 2 && !lineInternal::rightTurn(P[I[ml-2]], P[I[ml-1]], P[i]))
	popLeft();
      pushLeft(i);
    } else {
      while (1) {
	if (sizeRight() >= 2 && !lineInternal::leftTurn(P[I[mr+2]], P[I[mr+1]], P[i]))
	  popRight();
	else if (sizeRight() == 1 && !lineInternal::leftTurn(P[I[0]], P[I[mr+1]], P[i]))
	  popRight();
	else break;
      }
      pushRight(i);
    }
  }
  while (sizeLeft() >= 2 && !lineInternal::rightTurn(P[I[ml-2]], P[I[ml-1]], P[I[mr+1]]))
    popLeft();
#ifndef SILENT
  cout << " scan-time = " << t1.next() << endl;
#endif

  auto AR = I+n-sizeRight();
  for(intT i=0; i<sizeRight(); ++i) {
    I[sizeLeft()+i] = AR[i];
  }
#ifndef SILENT
  cout << " move-time = " << t1.stop() << endl;
#endif

  return sizeLeft()+sizeRight();
}

intT lineSweepExternalSerial(point2d* P, intT* Idx, intT n, intT* I) {
  static const bool verbose = false;

  if (!I || !Idx) abort();

  auto yCmp = [&](intT i, intT j) {
    point2d a = P[i];
    point2d b = P[j];
    return a.y() < b.y();
  };
  quickSortSerial(&Idx[0], n, yCmp);

  intT ml = 0;
  intT mr = n-1;

  auto pushLeft = [&](intT i) {I[ml++] = i;};
  auto popLeft = [&]() {ml--;};
  auto pushRight = [&](intT i) {I[mr--] = i;};
  auto popRight = [&]() {mr++;};
  auto sizeLeft = [&]{return ml;};
  auto sizeRight = [&]{return n-1-mr;};
  pushLeft(Idx[0]);

  //scan from bottom up
  for (intT i=1; i<n; ++i) {
    if (triArea(P[Idx[0]], P[Idx[n-1]], P[Idx[i]]) > numericKnob) {
      while (sizeLeft() >= 2 && !lineInternal::rightTurn(P[I[ml-2]], P[I[ml-1]], P[Idx[i]]))
	popLeft();
      pushLeft(Idx[i]);
    } else {
      while (1) {
	if (sizeRight() >= 2 && !lineInternal::leftTurn(P[I[mr+2]], P[I[mr+1]], P[Idx[i]]))
	  popRight();
	else if (sizeRight() == 1 && !lineInternal::leftTurn(P[I[0]], P[I[mr+1]], P[Idx[i]]))
	  popRight();
	else break;
      }
      pushRight(Idx[i]);
    }
  }
  while (sizeLeft() >= 2 && !lineInternal::rightTurn(P[I[ml-2]], P[I[ml-1]], P[I[mr+1]]))
    popLeft();

  auto AR = I+n-sizeRight();
  for(intT i=0; i<sizeRight(); ++i) {
    I[sizeLeft()+i] = AR[i];
  }

  return sizeLeft()+sizeRight();
}

#endif
