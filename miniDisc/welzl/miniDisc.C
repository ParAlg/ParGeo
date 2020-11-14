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

#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "miniDisc.h"
#include "geometry.h"
#include "check.h"

using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void miniDiscCaller(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef circle discT;
  static const bool serial = false;
  static const bool noRandom = true;
  static const bool moveToFront = true;
  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;

  if (n < 3) {
    cout << "smallest enclosing sphere needs at least 3 points" << endl;
    abort();}

  timing t0;t0.start();
  if(!noRandom) {
    cout << "permuting points" << endl;
    randPerm(P, n);
  }

  // auto p1 = point<3>(P[0].coordinate());
  // auto p2 = point<3>(P[1].coordinate());
  // auto p3 = point<3>(P[2].coordinate());
  // auto p4 = point<3>(P[3].coordinate());
  // cout << "p1 = " << p1 << endl;
  // cout << "p2 = " << p2 << endl;
  // cout << "p3 = " << p3 << endl;
  // cout << "p4 = " << p4 << endl;
  // cout << ">> sphere2" << endl;
  // sphere S = sphere(p1, p2);
  // cout << S.radius() << ", center " << S.center() << endl;
  // cout << ">> sphere3" << endl;
  // S = sphere(p1, p2, p3);
  // cout << S.radius() << ", center " << S.center() << endl;
  // cout << ">> sphere4" << endl;
  // S = sphere(p1, p2, p3, p4);
  // cout << S.radius() << ", center " << S.center() << endl;

  cout<< ">> grow ball" << endl;
  ball<dim> B = ball<dim>(P, 2);
  cout << B.radius() << ", center " << B.center() << endl;
  for (int i=2; i<dim+1; ++i) {
    B.grow(P[i]);
    cout << B.radius() << ", center " << B.center() << endl;
  }

  cout << "total-time = " << t0.stop() << endl;
  //cout << "disc = " << disc.center() << ", " << disc.radius() << endl;
  cout << endl;

  //check<2,discT>(&disc, P, n);
}

template<int dim>
void miniDisc(point<dim>* P, intT n) {
  return miniDiscCaller(P, n);
}

template void miniDisc<2>(point<2>*, intT);
template void miniDisc<3>(point<3>*, intT);
template void miniDisc<4>(point<4>*, intT);
template void miniDisc<5>(point<5>*, intT);
template void miniDisc<6>(point<6>*, intT);
template void miniDisc<7>(point<7>*, intT);
template void miniDisc<8>(point<8>*, intT);
template void miniDisc<9>(point<9>*, intT);
