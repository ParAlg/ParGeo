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
#include "bench.h"
#include "geometry.h"
using namespace std;

template<int dim>
sphere<dim> smallestEnclosingDisc2DSerial(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  auto circle = circleT(P[0].average(P[1]), P[0].pointDist(P[1]));
  intT i = 2;

  while (i < n) {
    bool conflict = false;
    for (; i<n; ++i) {
      if (!circle.contain(P[i])) {
        conflict = true;
        break;}
    }

    if (conflict) {
      cout << "conflict = " << i << endl;
      auto p0 = P[0];
      auto pi = P[i];//i now points to conflicting point
      circle = circleT(p0.average(pi), p0.pointDist(pi));
      for (intT j=1; j<i; ++j) {
        if (!circle.contain(P[j])) {
          circle = circleT(p0, pi, P[j]);
        }}
    } else i++;
  }
  return circle;
}

template<int dim>
sphere<dim> smallestEnclosingDisc2DParallel(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  auto flag = newA(intT, n+1);

  auto circle = circleT(P[0].average(P[1]), P[0].pointDist(P[1]));
  intT i = 2;

  while (i < n) {
    intT prefix = min(i+i, n);

    par_for (intT ii=i; ii<prefix; ++ii) {
      if (!circle.contain(P[ii])) flag[ii-i] = 1;
      else flag[ii-i] = 0;
    }
    intT numBad = sequence::prefixSum(flag, 0, prefix-i);
    flag[prefix-i] = numBad;

    if (numBad > 0) {
      intT conflict = -1;
      par_for(intT ii=0; ii<prefix-i; ++ii) {
        if (flag[ii]==0 && flag[ii]!=flag[ii+1]) conflict = ii;}
      i += conflict;//i now points to conflicting point
      cout << "conflict = " << i << endl;

      auto p0 = P[0];
      auto pi = P[i];
      circle = circleT(p0.average(pi), p0.pointDist(pi));
      for (intT j=1; j<i; ++j) {
        if (!circle.contain(P[j])) {
          circle = circleT(p0, pi, P[j]);
        }}
    } else {
      i = prefix;
    }
  }

  free(flag);
  return circle;
}

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  static const bool serial = false;
  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;
  if (dim > 2) {
    cout << "dim > 2 is not supported yet, abort" << endl;
    abort();}
  if (n < 2) abort();

  timing t0;t0.start();
  sphere<dim> circle = sphere<dim>();
  if(serial) circle = smallestEnclosingDisc2DSerial(P, n);
  else circle = smallestEnclosingDisc2DParallel(P, n);

  cout << "total-time = " << t0.stop() << endl;
  cout << "circle = " << circle.center << ", " << circle.radius << endl;

  //code for verifying correctness
  // for (intT i=0; i<n; ++i) {
  //   if(!circle.contain(P[0])) {
  //     cout << "outside point = " << P[0] << endl;
  //     cout << "dist = " << P[0].pointDist(circle.center) << endl;
  //     abort();
  //   }
  // }
  // cout << "correctness verified" << endl;
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
