// Copyright (c) 2020 Yiqiu Wang
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

#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "bench.h"
#include "geometry.h"
#include "miniDisc2d.h"
#include "miniDisc2dMTF.h"
//#include "miniDisc2dMTF2.h"
#include "miniDisc2dPrefix.h"
#include "miniDisc3d.h"
#include "check.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

void bench2D(point<2>* P, intT n) {
  typedef point<2> pointT;
  typedef circle discT;
  static const bool serial = false;
  static const bool noRandom = true;
  static const bool moveToFront = true;
  cout << "smallest enclosing disc, " << n << ", dim 2 points" << endl;

  if (n < 3) {
    cout << "smallest enclosing sphere needs at least 3 points" << endl;
    abort();}

  timing t0;t0.start();
  if(!noRandom) {
    cout << "permuting points" << endl;
    randPerm(P, n);
  }

  discT disc = discT();
  if(serial) {
    if (moveToFront) disc = MTF::miniDisc2DSerial(P, n);
    else disc = miniDisc2DSerial(P, n);
  } else {
    disc = miniDisc2DParallel(P, n);
    //disc = miniDisc2DPrefix(P, n);
  }
  cout << "total-time = " << t0.stop() << endl;
  cout << "disc = " << disc.center() << ", " << disc.radius() << endl;
  cout << endl;

  check<2,discT>(&disc, P, n);
}

void bench3D(point<3>* P, intT n) {
  typedef point<3> pointT;
  typedef sphere discT;
  static const bool serial = false;
  static const bool noRandom = true;
  cout << "smallest enclosing disc, " << n << ", dim 2 points" << endl;

  if (n < 4) {
    cout << "smallest enclosing sphere needs at least 4 points" << endl;
    abort();}

  timing t0;t0.start();
  if(!noRandom) {
    cout << "permuting points" << endl;
    randPerm(P, n);
  }

  discT disc = discT();
  if(serial) disc = miniDisc3DSerial(P, n);
  else disc = miniDisc3DParallel(P, n);
  cout << "total-time = " << t0.stop() << endl;
  cout << "disc = " << disc.center() << ", " << disc.radius() << endl;
  check<3,discT>(&disc, P, n);
}
