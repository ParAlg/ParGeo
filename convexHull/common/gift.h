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

#ifndef GIFT_H
#define GIFT_H

#include "geometry.h"
#include "hull.h"

namespace giftInternal {
  inline floatT angle(point2d& a, point2d& b, point2d& c) {
    point2d ab = b-a;
    point2d bc = c-b;
    return acos(ab.dot(bc) / (ab.length()*bc.length()));//this is SUPER expensive, todo
  };
}

intT giftWrapSerial(point2d* P, intT n, intT* I) {
  static const bool verbose = false;

  if (!I) abort();

  intT m=0;

  auto findLeft = [&](intT i) {return P[i].x();};
  intT si = sequence::minIndexSerial<floatT>(0, n, findLeft);
  point2d s = P[si];
  I[m++] = si;

  auto sp = point2d(s.x(), s.y()-1);
  auto findFirst = [&](intT j) {return giftInternal::angle(sp, s, P[j]);};
  intT myI = sequence::minIndexSerial<floatT>(0, n, findFirst);
  I[m++] = myI;
  if (verbose) {
    cout << "initial hull = " << P[I[0]] << " " << P[I[1]] << endl;
  }

  while (1) {
    auto a = P[I[m-2]];
    auto b = P[I[m-1]];
    auto findNext = [&](intT j) {return giftInternal::angle(a, b, P[j]);};
    myI = sequence::minIndex<floatT>(0, n, findNext);
    if (myI == I[0]) {
      break;
      if(verbose) cout << "hull closed" << endl;
    }
    I[m++] = myI;
  }
  return m;
}

intT giftWrapParallel(point2d* P, intT n, intT* I) {
  static const bool verbose = false;

  if (!I) abort();

  intT m=0;

  auto findLeft = [&](intT i) {return P[i].x();};
  intT si = sequence::minIndex<floatT>(0, n, findLeft);
  point2d s = P[si];
  I[m++] = si;

  auto sp = point2d(s.x(), s.y()-1);
  auto findFirst = [&](intT j) {return giftInternal::angle(sp, s, P[j]);};
  intT myI = sequence::minIndex<floatT>(0, n, findFirst);
  I[m++] = myI;

  if (verbose)
    cout << "initial hull = " << P[I[0]] << " " << P[I[1]] << endl;

  while (1) {
    auto a = P[I[m-2]];
    auto b = P[I[m-1]];
    auto findNext = [&](intT j) {return giftInternal::angle(a, b, P[j]);};
    myI = sequence::minIndex<floatT>(0, n, findNext);
    if (myI == I[0]) {
      break;
      if(verbose) cout << "hull closed" << endl;
    }
    I[m++] = myI;
  }
  return m;
}

#endif
