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

#ifndef MINI_DISC_2D
#define MINI_DISC_2D

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"
#include "check.h"
#include "prefix.h"

template<int dim>
sphere<dim> miniDisc2DSerial(point<dim>* P, intT n, point<dim> pi, point<dim> pj) {
  typedef sphere<dim> circleT;
  auto circle = circleT(pi, pj);

  for (intT i=0; i<n; ++i) {
    if (!circle.contain(P[i])) {
      circle = circleT(pi, pj, P[i]);
    }}
  return circle;
}

template<int dim>
sphere<dim> miniDisc2DSerial(point<dim>* P, intT n, point<dim> pi) {
  typedef sphere<dim> circleT;

  auto circle = circleT(P[0], pi);
  for (intT j=1; j<n; ++j) {
    if (!circle.contain(P[j])) {
      circle = miniDisc2DSerial(P, j, pi, P[j]);}
  }
  return circle;
}

template<int dim>
sphere<dim> miniDisc2DSerial(point<dim>* P, intT n) {
  typedef sphere<dim> circleT;

  auto circle = circleT(P[0], P[1]);
  for (intT i=2; i<n; ++i) {
    if (!circle.contain(P[i])) {
      cout << "ci = " << i << endl;
      circle = miniDisc2DSerial(P, i, P[i]);}
  }
  return circle;
}

//take ``vertical'' line (i,j) and get ``left/right'' most circle center
template<int dim>
sphere<dim> miniDisc2DParallel2(point<dim>* P, intT n, point<dim> pi, point<dim> pj, intT* flag=NULL) {
  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  auto circle = circleT(pi, pj);

  if (pi[1] == pj[1]) {
    abort();
  } else {
    bool freeFlag = false;
    if(!flag) {
      flag = newA(intT, n+1);//marks left/right to (pi,pj)
      freeFlag = true;}

    bool inputOk = true;
    par_for (intT jj=0; jj<n; ++jj) {
      if (!circle.contain(P[jj])) inputOk = false;
    }
    if (inputOk) return circle;

    pointT avg = pi.average(pj);
    floatT slope0 = (pi[1]-pj[1]) / (pi[0]-pj[0]);
    floatT offset0 = avg[1] - slope0*avg[0];

    auto centers = newA(pointT, n);
    par_for (intT jj=0; jj<n; ++jj) {
      auto c = circleT(pi, pj, P[jj]);

      //left/right
      floatT y = slope0*P[jj][0] + offset0;//y intersept with (pi,pj)
      if ((y < P[jj][1] && slope0 >= 0) || (y >= P[jj][1] && slope0 < 0)) {
        flag[jj] = 0;
      } else flag[jj] = 1;

      centers[jj] = c.center();
    }

    //left/right center extrema from left/right points
    struct getL {
      pointT* in; intT* flag;
      getL(pointT* inn, intT* flagg): in(inn), flag(flagg) {};
      floatT operator() (intT idx) {
        if (flag[idx] == 0) return in[idx][0];
        else return floatMax();
      }};
    struct getR {
      pointT* in; intT* flag;
      getR(pointT* inn, intT* flagg): in(inn), flag(flagg) {};
      floatT operator() (intT idx) {
        if (flag[idx] == 1) return in[idx][0];
        else return floatMin();
      }};

    intT iMin = sequence::minIndex<floatT, intT, getL>(0, n, getL(centers, flag));
    intT iMax = sequence::maxIndex<floatT, intT, getR>(0, n, getR(centers, flag));
    auto circle1 = circleT(pi, pj, P[iMin]);
    auto circle2 = circleT(pi, pj, P[iMax]);

    if(freeFlag) free(flag);
    free(centers);

    //among two, return the one that passes the check, and prefer smaller radius
    bool check1 = check<dim>(&circle1, P, n, false);
    if (check1 && circle1.radius() < circle2.radius()) return circle1;
    bool check2 = check<dim>(&circle2, P, n, false);
    if (check2 && circle2.radius() < circle1.radius()) return circle2;

    if (check1 && check2)
      return circle1.radius() < circle2.radius() ? circle1 : circle2;
    else if (check1) return circle1;
    else if (check2) return circle2;
    else {
      cout << "error, no valid circle, abort()" << endl;
      abort();}

  }
}

template<int dim>
sphere<dim> miniDisc2DParallel(point<dim>* P, intT n, point<dim> pi, point<dim> pj, intT* flag=NULL) {
  typedef sphere<dim> circleT;
  typedef point<dim> pointT;

  auto circle = circleT(pi, pj);
  auto process = [&](pointT p) {
                   if (!circle.contain(p)) return true;
                   else return false;
                 };
  auto cleanUp = [&](pointT* A, intT ci) {
                   circle = circleT(pi, pj, A[ci]);
                 };
  parallel_prefix(P, n, process, cleanUp, false, flag);
  return circle;
}

template<int dim>
sphere<dim> miniDisc2DParallel(point<dim>* P, intT n, point<dim> pi, intT* flag=NULL) {
  typedef sphere<dim> circleT;
  typedef point<dim> pointT;

  auto circle = circleT(P[0], pi);
  auto process = [&](pointT p) {
                   if (!circle.contain(p)) return true;
                   else return false;
                 };
  auto cleanUp = [&](pointT* A, intT ci) {
                   if (ci < 4000)
                     circle = miniDisc2DParallel(A, ci, pi, A[ci], flag);
                   else
                     circle = miniDisc2DParallel2(A, ci, pi, A[ci], flag);
                 };
  parallel_prefix(P, n, process, cleanUp, false, flag);
  return circle;
}

template<int dim>
sphere<dim> miniDisc2DParallel(point<dim>* P, intT n) {
  typedef sphere<dim> circleT;
  typedef point<dim> pointT;

  intT* flag = newA(intT, n+1);

  auto circle = circleT(P[0], P[1]);
  auto process = [&](pointT p) {
                   if (!circle.contain(p)) return true;
                   else return false;
                 };
  auto cleanUp = [&](pointT* A, intT ci) {
                   circle = miniDisc2DParallel(A, ci, A[ci], flag);
                 };
  parallel_prefix(P, n, process, cleanUp, true, flag);

  free(flag);
  return circle;
}

#endif
