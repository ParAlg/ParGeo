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

#ifndef MIN_DISC_2D
#define MIN_DISC_2D

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"
#include "check.h"

//update 2: search in P[0:j) a third point for circle(i,j) (where i < j)
// that makes the circle include all points <= i.

//dim=2
template<int dim>
sphere<dim> update2Serial(point<dim>* P, intT i, intT j, sphere<dim>& circle) {
  typedef sphere<dim> circleT;

  for (intT jj=0; jj<j; ++jj) {
    if (!circle.contain(P[jj])) {
      circle = circleT(P[i], P[j], P[jj]);
    }}
  return circle;
}

//s: vector through origin, offset: y offset, v: input
template<int dim>
point<dim> project(point<dim> s, floatT offset, point<dim> v) {
  v[1] = v[1]-offset;
  auto r = s.mult(v.dot(s)/s.dot(s));
  r[1] = r[1]+offset;
  return r;
}

//dim=2
//take ``vertical'' line (i,j) and get left/right most circle center
template<int dim>
sphere<dim> update2Parallel(point<dim>* P, intT i, intT j, sphere<dim> circle) {
  //return update2Serial(P, i, j, circle);
  if (j < 100) return update2Serial(P, i, j, circle);

  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  if (P[i][0] == P[j][0]) {
    return update2Serial(P, i, j, circle);
  } else if (P[i][1] == P[j][1]) {
    return update2Serial(P, i, j, circle);
  } else {
    cout << "points = " << j << endl;

    auto flag = newA(intT, 2*(j+1));
    auto flag2 = flag + j + 1;

    par_for (intT jj=0; jj<j; ++jj) {
      if (!circle.contain(P[jj])) flag[jj] = 1;
      else flag[jj] = 0;}
    intT numCirc = sequence::prefixSum(flag, 0, j);
    flag[j] = numCirc;
    cout << "circles = " << numCirc << endl;

    if (numCirc <= 0) return circle;

    pointT avg = P[i].average(P[j]);
    floatT slope0 = abs(P[i][1]-P[j][1]) / abs(P[i][0]-P[j][0]);
    floatT offset0 = avg[1] - slope0*avg[0];
    floatT slope = -1/slope0;//of bisector
    floatT offset = avg[1] - slope*avg[0];//bisector: y = slope*x + offset
    auto s = pointT(); s[0] = 1; s[1] = slope;
    cout << "P[i] = " << P[i] << endl;
    cout << "P[j] = " << P[j] << endl;
    // cout << "offset0 = " << offset0 << endl;
    // cout << "offset = " << offset << endl;

    auto circs = newA(circleT, numCirc);
    par_for (intT jj=0; jj<j; ++jj) {
      if(flag[jj] != flag[jj+1]) {//form by a point not in circ(P[i], P[j])
        auto c = circleT(P[i], P[j], P[flag[jj]]);
        circs[flag[jj]] = c;
        floatT y = slope0*c.center()[0] + offset0;//y intersept with (P[i],P[j])
        if ((y < c[1] && slope0 >= 0) || (y >= c[1] && slope0 < 0)) {
          //center is left of (P[i],P[j])
          flag2[flag[jj]] = 0;
        } else {
          //center is right of (P[i],P[j])
          flag2[flag[jj]] = 1;
        }
      }
    }

    //projection
    auto proj = newA(pointT, numCirc);
    par_for (intT jj=0; jj<j; ++jj) {
      if(flag[jj] != flag[jj+1]) {
        proj[flag[jj]] = project(s, offset, circs[flag[jj]].center());
      }}

    //left/right extrema, find min/max along x
    struct getL {//max
      pointT* in; intT* flag;
      getL(pointT* inn, intT* flagg): in(inn), flag(flagg) {};
      floatT operator() (intT idx) {
        if (flag[idx] == 0) return in[idx][0];
        else return floatMin();
      }
    };
    struct getR {//min
      pointT* in; intT* flag;
      getR(pointT* inn, intT* flagg): in(inn), flag(flagg) {};
      floatT operator() (intT idx) {
        if (flag[idx] == 1) return in[idx][0];
        else return floatMax();
      }
    };

    intT iMin = sequence::minIndex<floatT, intT, getL>(0, numCirc, getL(proj, flag2));
    intT iMax = sequence::maxIndex<floatT, intT, getR>(0, numCirc, getR(proj, flag2));
    auto circle1 = circs[iMin];
    auto circle2 = circs[iMax];

    cout << "circ1 = " << circle1.center() << ", " << circle1.radius() << endl;
    cout << "circ2 = " << circle2.center() << ", " << circle2.radius() << endl;
    cout << "check 1 " << check<dim>(&circle1, P, j, false) << endl;

    cout << "check 2 " << check<dim>(&circle2, P, j, false) << endl;
    cout << "check " << check<dim>(&circle, P, j, false) << endl;
    abort();//todo, debug, implementation is problematic

    check<dim>(&circle2, P, j);
    free(flag);
    free(circs);
    free(proj);
    cout << endl;

    //take min radius among three circles
    if (circle.radius() < circle1.radius()) {
      return circle.radius() < circle2.radius() ? circle : circle2;
    } else {
      return circle1.radius() < circle2.radius() ? circle1 : circle2;
    }
  }
}

template<int dim>
sphere<dim> minDisc2DSerial(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  auto circle = circleT(P[0], P[1]);
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

      circle = circleT(P[0], P[i]);

      for (intT j=1; j<i; ++j) {//update 1
        if (!circle.contain(P[j])) {
          circle = circleT(P[i], P[j]);
          circle = update2Serial(P, i, j, circle);//update 2
        }}

    } else i++;
  }
  return circle;
}

template<int dim>
sphere<dim> minDisc2DParallel(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef sphere<dim> circleT;

  auto flag = newA(intT, n+1);

  auto circle = circleT(P[0], P[1]);
  intT i = 2;

  while (i < n) {
    intT prefix = min(i+i, n);

    bool parInsertion = false;
    bool conflict = false;

    //insertion
    if (prefix - i < 2000) {
      //serial insertion
      for (; i<prefix; ++i) {
        if (!circle.contain(P[i])) {
          conflict = true;
          break;}
      }
    } else {
      //parallel insertion
      parInsertion = true;
      par_for (intT ii=i; ii<prefix; ++ii) {
        if (!circle.contain(P[ii])) flag[ii-i] = 1;
        else flag[ii-i] = 0;
      }
      intT numBad = sequence::prefixSum(flag, 0, prefix-i);
      flag[prefix-i] = numBad;
      if (numBad > 0) conflict = true;
    }

    if (conflict) {
      if (parInsertion) {
        intT conflict = -1;
        par_for(intT ii=0; ii<prefix-i; ++ii) {
          if (flag[ii]==0 && flag[ii]!=flag[ii+1]) conflict = ii;}
        i += conflict;//i now points to conflicting point
      }
      cout << "conflict = " << i << endl;

      circle = circleT(P[0], P[i]);

      //update 1
      if (i < 2000) {
        //serial update 1
        for (intT j=1; j<i; ++j) {
          if (!circle.contain(P[j])) {
            circle = circleT(P[i], P[j]);
            circle = update2Serial(P, i, j, circle);
          }
        }
      } else {
        //parallel update 1
        intT j = 1;
        while (j<i) {
          intT prefix = min(j+j, i);
          intT conflict = -1;
          par_for (intT jj=j; jj<prefix; ++jj) {
            if (!circle.contain(P[jj])) flag[jj-j] = 1;
            else flag[jj-j] = 0;
          }
          intT numBad = sequence::prefixSum(flag, 0, prefix-j);
          flag[prefix-j] = numBad;

          if (numBad > 0) {
            par_for(intT jj=0; jj<prefix-j; ++jj) {
              if (flag[jj]==0 && flag[jj]!=flag[jj+1]) conflict = jj;}
            j += conflict;

            circle = circleT(P[i], P[j]);
            circle = update2Serial(P, i, j, circle);
          } else {
            j = prefix;
          }
        }
      }

    } else {
      i = prefix;
    }
  }

  free(flag);
  return circle;
}
#endif
