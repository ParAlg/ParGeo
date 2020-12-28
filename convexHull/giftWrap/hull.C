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
#include "pbbs/utils.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "hull.h"
#include "serialQuick.h"
using namespace std;

_seq<intT> hull(point2d* P, intT n) {
  static const bool verbose = false;
  auto angle = [&](point2d& a, point2d& b, point2d& c) {
		 point2d ab = b-a;
		 point2d bc = c-b;
		 return acos(ab.dot(bc) / (ab.length()*bc.length()));
	       };

  timing t; t.start();
  intT* I = newA(intT, n);
  for (intT i=0; i < n; i++) I[i] = i;
  intT m=0;

  //find left-most point
  point2d s = P[0];
  intT si = 0;
  for (intT i=0; i < n; i++) {
    if (P[i].x() < s.x()) {
      s = P[i];
      si = i;
    }
  }
  I[m++] = si;

  auto sp = point2d(s.x(), s.y()-1);
  floatT myMin = floatMax();
  intT myI = -1;
  for(intT j=0; j<n; ++j) {
    floatT ag = angle(sp, s, P[j]);
    if (ag < myMin) {
      myMin = ag;
      myI = j;}
  }
  I[m++] = myI;
  if (verbose) {
    cout << "initial hull = " << P[I[0]] << " " << P[I[1]] << endl;
    cout << "angle = " << myMin << ", size = " << m << endl;
  }

  while (1) {
    floatT myMin = floatMax();
    intT myI = -1;
    auto a = P[I[m-2]];
    auto b = P[I[m-1]];
    for(intT j=0; j<n; ++j) {
      floatT ag = angle(a, b, P[j]);
      if (ag < myMin) {
	myMin = ag;
	myI = j;}
    }
    if (myI == I[0]) {
      break;
      if(verbose) cout << "hull closed" << endl;
    }
    I[m++] = myI;
    if(verbose) cout << "angle = " << myMin << ", size = " << m << endl;
  }

  cout << "hull-time = " << t.stop() << endl;

  check(P, n, I, m);
  return _seq<intT>(I, m);
}
