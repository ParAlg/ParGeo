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

#ifndef MINI_DISC_2D_MTF2_H
#define MINI_DISC_2D_MTF2_H

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"
#include "check.h"
#include "prefix.h"

#include <iostream>
#include <fstream>
#include <string>

static const bool mtf = true;
intT confCount = 0;

namespace MTF {
  circle miniDisc2DSerial(point<2>* P, point<2>* endPtr, point<2> pi, point<2> pj) {
    typedef circle circleT;

    auto circle = circleT(pi, pj);

    auto ptr = P;
    while (ptr < endPtr) {
      auto p = *(ptr);
      if (p.isEmpty()) {
        ptr++;
        continue;}
      if (!circle.contain(p)) {
        confCount ++;
        circle = circleT(pi, pj, p);
        if (mtf) {
          //?
        }
      }
      ptr++;
    }
    return circle;
  }

  circle miniDisc2DSerial(point<2>* P, point<2>* endPtr, point<2> pi) {
    typedef circle circleT;
    typedef point<2> pointT;

    //find first non-empty point
    intT p0 = 0;
    while (P[p0].isEmpty()) p0++;
    auto circle = circleT(P[p0], pi);

    auto ptr = P;
    while (ptr < endPtr) {
      auto p = *(ptr);
      if (p.isEmpty()) {
        ptr++;
        continue;}
      if (!circle.contain(p)) {
        confCount ++;
        circle = miniDisc2DSerial(P, ptr, pi, p);
        if (mtf) {
          // ?
        }
      }
      ptr++;
    }

    return circle;
  }

  circle miniDisc2DSerial(point<2>* PP, intT n) {
    typedef circle circleT;
    typedef point<2> pointT;

    pointT* P;
    pointT* PP2;
    if (mtf) {
      PP2 = newA(pointT, n*2);
      for(intT i=0; i<n; ++i) PP2[n+i] = PP[i];
      P = PP2+n;
    } else {
      P = PP;
    }

    auto circle = circleT(P[0], P[1]);

    auto ptr = P;
    auto endPtr = P+n;
    while (ptr < endPtr) {
      auto p = *(ptr);
      if (p.isEmpty()) {
        ptr++;
        continue;}
      if (!circle.contain(p)) {
        confCount ++;
        circle = miniDisc2DSerial(P, ptr, p);
        if (mtf) {
          ptr[0] = pointT();
          P--;
          P[0] = p;
        }
      }
      ptr++;
    }

    cout << "#conflicts = " << confCount << endl;
    if (mtf) free(PP2);
    return circle;
  }

}

#endif
