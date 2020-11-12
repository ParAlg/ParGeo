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

#ifndef MINI_DISC_2D_MTF_H
#define MINI_DISC_2D_MTF_H

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"
#include "check.h"
#include "prefix.h"

#include <iostream>
#include <fstream>
#include <string>

static const bool takeSnapshot = false;
static const bool mtf = true;
intT ii = 0;
intT confCount = 0;

void snapshot(const char* name, point<2>* P, intT n, intT ci, circle c, point<2> pi, intT idxi, point<2> pj, intT idxj) {
  cout << "snapshot-lvl3 " << ii-1 << ": n = " << n << ", pi = " << idxi << ", pj = " << idxj << endl;
  std::ofstream myfile;
  myfile.open(name);
  myfile << c.center() << "";
  myfile << c.radius() << endl;
  myfile << ci << endl;
  for (intT i=0; i<n; ++i) myfile << P[i] << "" << i << endl;
  myfile << pi << "" << idxi << endl;
  myfile << pj << "" << idxj << endl;
  myfile.close();
}

void snapshot(const char* name, point<2>* P, intT n, intT ci, circle c, point<2> pi, intT idxi) {
  cout << "snapshot-lvl2 " << ii-1 << ": n = " << n << ", pi = " << idxi << endl;
  std::ofstream myfile;
  myfile.open(name);
  myfile << c.center() << "";
  myfile << c.radius() << endl;
  myfile << ci << endl;
  for (intT i=0; i<n; ++i) myfile << P[i] << "" << i << endl;
  myfile << pi << "" << idxi << endl;
  myfile.close();
}

void snapshot(const char* name, point<2>* P, intT n, intT ci, circle c) {
  cout << "snapshot-lvl1 " << ii-1 << ": n = " << n << endl;
  std::ofstream myfile;
  myfile.open(name);
  myfile << c.center() << "";
  myfile << c.radius() << endl;
  myfile << ci << endl;
  for (intT i=0; i<n; ++i) myfile << P[i] << "" << i << endl;
  myfile.close();
}

namespace Heuristic {

  void moveToFront(point<2>* P, intT j, point<2>* space) {
    // cout << "swap 0 and " << j << endl;
    // cout << "before = ";for(intT i=0; i<10; ++i) cout << P[i] << " ";cout << endl;
    par_for(intT i=0; i<=j; ++i) space[i] = P[i];
    P[0] = space[j];
    par_for(intT i=1; i<=j; ++i) P[i] = space[i-1];
    // cout << "after = ";for(intT i=0; i<10; ++i) cout << P[i] << " ";cout << endl;
  }

  circle miniDisc2DSerial(point<2>* P, intT n, point<2> pi, intT idxi, point<2> pj, intT idxj, point<2>* PP) {
    typedef circle circleT;
    auto circle = circleT(pi, pj);

    for (intT i=0; i<n; ++i) {
      if (takeSnapshot) {
        auto fname = "./snapshot/" + std::to_string(ii++) + ".txt";
        snapshot(fname.c_str(), P, n, i, circle, pi, idxi, pj, idxj);}
      if (!circle.contain(P[i])) {
        confCount ++;
        circle = circleT(pi, pj, P[i]);
        if (mtf) {
          swap(P[0], P[i]);
          //moveToFront(P, i, PP);
        }
      }
    }
    return circle;
  }

  circle miniDisc2DSerial(point<2>* P, intT n, point<2> pi, intT idxi, point<2>* PP) {
    typedef circle circleT;

    auto circle = circleT(P[0], pi);
    for (intT j=0; j<n; ++j) {
      if (takeSnapshot) {
        auto fname = "./snapshot/" + std::to_string(ii++) + ".txt";
        snapshot(fname.c_str(), P, n, j, circle, pi, idxi);}
      if (!circle.contain(P[j])) {
        confCount ++;
        circle = miniDisc2DSerial(P, j, pi, idxi, P[j], j, PP);
        if (mtf) {
          swap(P[1], P[j]);
          //moveToFront(P, j, PP);
        }
      }
    }
    return circle;
  }

  circle miniDisc2DSerialMTF(point<2>* P, intT n) {
    typedef circle circleT;
    typedef point<2> pointT;

    pointT* PP = newA(pointT, n);

    auto circle = circleT(P[0], P[1]);

    auto findPivot = [&] (intT s)
                     {
                       floatT rSqr = circle.radius() * circle.radius();
                       floatT dMax = 0;
                       intT bestI = -1;
                       for (intT ii=s; ii<n; ++ii) {
                         floatT tmp = P[ii].distSqr(circle.center());
                         if (tmp - rSqr > dMax) { // ||p-c||^2 - r^2
                           bestI = ii;
                           dMax = tmp - rSqr;
                         }
                       }
                       return bestI;
                     };

    for (intT i=0; i<n; ++i) {
      if (takeSnapshot) {
        auto fname = "./snapshot/" + std::to_string(ii++) + ".txt";
        snapshot(fname.c_str(), P, n, i, circle);}

      if (!circle.contain(P[i])) {
        //pivoting
        intT ii = findPivot(i);
        swap(P[ii],P[i]);

        confCount ++;
        cout << "ci = " << i << ", " << P[i] << endl;
        circle = miniDisc2DSerial(P, i, P[i], i, PP);
        if (mtf) {
          swap(P[2], P[i]);
          //moveToFront(P, i, PP);
        }
      }
    }

    cout << "#conflicts = " << confCount << endl;
    free(PP);
    return circle;
  }

  circle miniDisc2DSerialPivot(point<2>* P, intT n) {
    typedef circle circleT;
    typedef point<2> pointT;
    pointT* PP = newA(pointT, n);

    auto circle = circleT(P[0], P[1]);

    auto findPivot = [&] ()
                     {
                       floatT rSqr = circle.radius() * circle.radius();
                       floatT dMax = 0;
                       intT bestI = -1;
                       for (intT ii=0; ii<n; ++ii) {
                         floatT tmp = P[ii].distSqr(circle.center());
                         if (tmp - rSqr > dMax) { // ||p-c||^2 - r^2
                           bestI = ii;
                           dMax = tmp - rSqr;
                         }
                       }
                       return bestI;
                     };

    while(1) {
      intT i = findPivot();
      if (i < 0) break;

      if (takeSnapshot) {
        auto fname = "./snapshot/" + std::to_string(ii++) + ".txt";
        snapshot(fname.c_str(), P, n, i, circle);}

      if (!circle.contain(P[i])) {
        confCount ++;
        cout << "ci = " << i << ", " << P[i] << endl;
        circle = miniDisc2DSerial(P, i, P[i], i, PP);
        if (mtf) swap(P[2], P[i]);
      }
    }

    cout << "#conflicts = " << confCount << endl;
    free(PP);
    return circle;
  }

}

#endif
