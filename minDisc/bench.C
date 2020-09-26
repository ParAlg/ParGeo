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
#include "minDisc2d.h"
#include "minDisc3d.h"
#include "check.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  static const bool serial = false;
  static const bool noRandom = true;
  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;
  if (dim > 3) {
    cout << "smallest enclosing sphere only supported for dim <= 3" << endl;
    abort();}

  if (n < 4) {
    cout << "smallest enclosing sphere needs at least 4 points" << endl;
    abort();}

  timing t0;t0.start();
  if(!noRandom) {
    //permutation
    cout << "permuting points" << endl;
    randPerm(P, n);
  }

  if (dim == 2) {
    sphere<dim> disc = sphere<dim>();
    if(serial) disc = minDisc2DSerial(P, n);
    else disc = minDisc2DParallel(P, n);
    cout << "total-time = " << t0.stop() << endl;
    cout << "disc = " << disc.center() << ", " << disc.radius() << endl;
    check(&disc, P, n);
  } else if (dim == 3) {
    sphere<dim> disc = sphere<dim>();
    if(serial) disc = minDisc3DSerial(P, n);
    else disc = minDisc3DParallel(P, n);
    cout << "total-time = " << t0.stop() << endl;
    cout << "disc = " << disc.center() << ", " << disc.radius() << endl;
    check(&disc, P, n);
  } else {
    cout << "smallest enclosing sphere only supported for dim <= 3" << endl;
    abort();
  }

}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
