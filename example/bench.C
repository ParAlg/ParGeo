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
#include "pbbs/sequence.h"
#include "pbbs/randPerm.h"
#include "bench.h"
#include "geometry.h"
using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  cout << "testing points, " << n << ", dim " << dim << " points" << endl;
  if (dim != 3) {
    cout << "this test needs a dim==3 dataset, abort" << endl;
    abort();}
  if (n < 4) abort();

  auto p0 = point3d(P[0].coordinate());
  auto p1 = point3d(P[1].coordinate());
  auto p2 = point3d(P[2].coordinate());
  auto v0 = vect3d(P[0].coordinate());
  auto v1 = vect3d(P[1].coordinate());
  auto v2 = vect3d(P[2].coordinate());

  cout << "p0 = " << p0 << endl;
  cout << "p1 = " << p1 << endl;
  cout << "p2 = " << p2 << endl;
  cout << "v0 = " << v0 << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;

  cout << "p1 + p2 = " << p1+p2 << endl;
  cout << "p1 - p2 = " << p1-p2 << endl;
  cout << "p1.minCoords(p2) = " << p1.minCoords(p2) << endl;
  cout << "p1.maxCoords(p2) = " << p1.maxCoords(p2) << endl;
  cout << "p0 + p1 - p2 = " << p0+p1-p2 << endl;

  cout << "v1 + v2 = " << v1+v2 << endl;
  cout << "v1 + p2 = " << v1+p2 << endl;
  cout << "p2 + v1 = " << p2+v1 << endl;
  cout << "v1 - v2 = " << v1-v2 << endl;
  cout << "v1.dot(v2) = " << v1.dot(v2) << endl;
  cout << "v1.cross(v2) = " << v1.cross(v2) << endl;
  cout << "v1 / 10 = " << v1/10 << endl;
  cout << "v0 + v1 - v2 = " << v0+v1-v2 << endl;

  cout << endl << "<<<<<<<<<<<<<<<" << endl;
  auto pp0 = point2d(P[0].coordinate());
  auto pp1 = point2d(P[1].coordinate());
  auto pp2 = point2d(P[2].coordinate());
  auto pp3 = point2d(P[3].coordinate());
  auto vv0 = vect2d(P[0].coordinate());
  cout << "pp0 = " << pp0 << endl;
  cout << "pp1 = " << pp1 << endl;
  cout << "pp2 = " << pp2 << endl;
  cout << "pp3 = " << pp3 << endl;
  cout << "triArea(pp0, pp1, pp2) = " << triArea(pp0, pp1, pp2) << endl;
  cout << "triAreaNormalized(pp0, pp1, pp2) = " << triAreaNormalized(pp0, pp1, pp2) << endl;
  cout << "counterClockwise(pp0, pp1, pp2) = " << counterClockwise(pp0, pp1, pp2) << endl;
  cout << "onParabola(vv0) = " << onParabola(vv0) << endl;
  cout << "inCircle(pp0, pp1, pp2, pp3) = " << inCircle(pp0, pp1, pp2, pp3) << endl;
  cout << "inCircleNormalized(pp0, pp1, pp2, pp3) = " << inCircleNormalized(pp0, pp1, pp2, pp3) << endl;
  cout << "angle(pp0, pp1, pp2) = " << angle(pp0, pp1, pp2) << endl;
  cout << "minAngle(pp0, pp1, pp2) = " << minAngle(pp0, pp1, pp2) << endl;
  cout << "minAngleCheck(pp0, pp1, pp2, 10) = " << minAngleCheck(pp0, pp1, pp2, 10) << endl;
  cout << "triangleCircumcenter(pp0, pp1, pp2) = " << triangleCircumcenter(pp0, pp1, pp2) << endl;
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
