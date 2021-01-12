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

#ifndef GRAHAM_H
#define GRAHAM_H

#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "hull.h"

_seq<intT> grahamScanSerial(point2d* P, intT n, intT* I=NULL) {
  static const bool verbose = false;
  auto angle = [&](point2d& a, point2d& b, point2d& c) {
    point2d ab = b-a;
    point2d bc = c-b;
    return acos(ab.dot(bc) / (ab.length()*bc.length()));
  };
  static const intT left = 0;
  static const intT right = 1;
  static const intT straight = 2;
  auto turn = [&](point2d& a, point2d& b, point2d& c) {
    auto cross = (b.x()-a.x())*(c.y()-a.y()) - (b.y()-a.y())*(c.x()-a.x());
    if (cross > 0) return left;
    else if (cross < 0) return right;
    else return straight;
  };

  if (!I) {
    I = newA(intT, n);
  }
  intT m=0;
#ifndef SILENT
  timing t; t.start();
#endif
  auto findLeft = [&](intT i) {return P[i].x();};
  intT si = sequence::minIndexSerial<floatT>(0, n, findLeft);
  point2d s = P[si];
  //point2d sbelow = s;
  //sbelow[1] = s.y()-1;
#ifndef SILENT
  cout << "init-time = " << t.next() << endl;
#endif
  auto sp = point2d(s.x(), s.y()-1);
  auto angleLess = [&](point2d a, point2d b) {
    //return angle(sbelow, s, a) < angle(sbelow, s, b);
    intT myTurn = turn(s, a, b);
    if (myTurn == right) return true;
    else if (myTurn == left) return false;
    else return a.x() < b.x();
  };
  swap(P[0], P[si]);
  sampleSort(&P[1], n-1, angleLess);
  //sort(&P[1], &P[n-1], angleLess);

  I[m++] = 0;
  I[m++] = 1;
#ifndef SILENT
  cout << "sort-time = " << t.next() << endl;
#endif
  auto push = [&](intT i) {I[m++] = i;};
  auto pop = [&]() {m--;};
  auto isEmpty = [&]() {return m==0;};

  for (intT i=2; i<n; ++i) {
    if (turn(P[I[m-2]], P[I[m-1]], P[i])!=right) {
      pop();
      while (turn(P[I[m-2]], P[I[m-1]], P[i])!=right) pop();
    }
    push(i);
  }
#ifndef SILENT
  cout << "scan-time = " << t.stop() << endl;
#endif
  return _seq<intT>(I, m);
}

#endif
