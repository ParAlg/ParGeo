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

#ifndef MINI_DISC_3D
#define MINI_DISC_3D

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"

sphere miniDisc3DSerial(point<3>* P, intT n, point<3> pi, point<3> pj, point<3> pk) {
  typedef sphere discT;

  auto disc = discT(pi, pj, pk);
  for (intT l=0; l<n; ++l) {
    if (!disc.contain(P[l])) {
      disc = discT(pi, pj, pk, P[l]);}
  }
  return disc;
}

sphere miniDisc3DSerial(point<3>* P, intT n, point<3> pi, point<3> pj) {
  typedef sphere discT;

  auto disc = discT(pi, pj);
  for (intT k=0; k<n; ++k) {
    if (!disc.contain(P[k])) {
      disc = miniDisc3DSerial(P, k, pi, pj, P[k]);}
  }
  return disc;
}

sphere miniDisc3DSerial(point<3>* P, intT n, point<3> pi) {
  typedef sphere discT;

  auto disc = discT(P[0], pi);
  for (intT j=1; j<n; ++j) {
    if (!disc.contain(P[j])) {
      disc = miniDisc3DSerial(P, j, pi, P[j]);}
  }
  return disc;
}

sphere miniDisc3DSerial(point<3>* P, intT n) {
  typedef sphere discT;

  auto disc = discT(P[0], P[1]);
  for (intT i=2; i<n; ++i) {
    if (!disc.contain(P[i])) {
      cout << "ci = " << i << endl;
      disc = miniDisc3DSerial(P, i, P[i]);}
  }
  return disc;
}


sphere miniDisc3DParallel(point<3>* P, intT n) {
  miniDisc3DSerial(P, n);
}

#endif
