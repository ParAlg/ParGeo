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

#include <vector>
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "prefix.h"
#include "miniDisc.h"
#include "geometry.h"
#include "welzl.h"
#include "check.h"

using namespace std;

template<int dim>
ball<dim> miniDiscMtfSerial(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  B = support2Ball(P, support);

  if (B.size() == dim+1) {
    return B;
  }

  for (intT i=0; i<n; ++i) {

    if (!B.contain(P[i])) {
      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();
      support.push_back(P[i]);
      B = miniDiscMtfSerial(P, i, support, B);
      support.pop_back();

      if (i > dim-support.size()) {
        swap(P[dim-support.size()], P[i]);}
    }
  }

  return B;
}

template<int dim>
ball<dim> miniDiscMtf(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B, intT* flag=NULL) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  if (n < 2000) return miniDiscMtfSerial(P, n, support, B);

  B = support2Ball(P, support);
  if (B.size() == dim+1) return B;

  bool freeFlag = false;
  if (!flag) {
    freeFlag = true;
    flag = newA(intT, n+1);}

  auto process = [&](pointT p) {
                   if (!B.contain(p)) return true;
                   else return false;
                 };
  auto cleanUp = [&](pointT* A, intT i) {
                   if (support.size() == B.size()) B.grow(A[i]);
                   else B = ballT();
                   support.push_back(A[i]);
                   B = miniDiscMtf(A, i, support, B, flag);
                   support.pop_back();

                   if (i > dim-support.size()) {
                     swap(P[dim-support.size()], P[i]);}
                 };
  parallel_prefix(P, n, process, cleanUp, freeFlag, flag);

  if(freeFlag) free(flag);
  return B;
}

template<int dim>
intT findPivot(point<dim>* P, intT n, ball<dim> B, intT s) {
  floatT rSqr = B.radius() * B.radius();
  floatT dMax = 0;
  intT bestI = -1;
  for (intT ii=s; ii<n; ++ii) {
    floatT tmp = P[ii].distSqr(B.center());
    if (tmp - rSqr > dMax) { // ||p-c||^2 - r^2
      bestI = ii;
      dMax = tmp - rSqr;
    }
  }
  return bestI;
}

template<int dim>
ball<dim> miniDiscPivotSerial(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  B = support2Ball(P, support);

  if (B.size() == dim+1) {
    return B;
  }

  for (intT i=0; i<n; ++i) {
    if (!B.contain(P[i])) {
      intT ii = findPivot(P, n, B, i+1);
      if (ii > i) swap(P[ii], P[i]);

      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();
      support.push_back(P[i]);
      B = miniDiscMtfSerial(P, i, support, B);
      support.pop_back();

      if (i > dim-support.size()) {
        swap(P[dim-support.size()], P[i]);}
    }
  }

  return B;
}

template<int dim>
ball<dim> miniDiscPivot(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B, intT* flag=NULL) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  if (n < 2000) return miniDiscPivotSerial(P, n, support, B);

  B = support2Ball(P, support);
  if (B.size() == dim+1) return B;

  bool freeFlag = false;
  if (!flag) {
    freeFlag = true;
    flag = newA(intT, n+1);}

  auto process = [&](pointT p) {
                   if (!B.contain(p)) return true;
                   else return false;
                 };
  auto cleanUp = [&](pointT* A, intT i) {
                   intT ii = findPivot(A, i, B, i+1);
                   if (ii > i) swap(P[ii], P[i]);

                   if (support.size() == B.size()) B.grow(A[i]);
                   else B = ballT();
                   support.push_back(A[i]);
                   B = miniDiscMtf(A, i, support, B, flag);
                   support.pop_back();

                   if (i > dim-support.size()) {
                     swap(P[dim-support.size()], P[i]);}
                 };
  parallel_prefix(P, n, process, cleanUp, freeFlag, flag);

  if(freeFlag) free(flag);
  return B;
}

template<int dim>
void miniDisc(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef circle discT;
  typedef ball<dim> ballT;

  static const bool preprocess = false;
  static const bool serial = false;
  /*
    - 0: plain
    - 1: mtf
    - 2: pivot + mtf
   */
  static const int method = 2;

  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;

  timing t0;t0.start();
  if(preprocess) {
    randPerm(P, n);
    cout << "preprocess-time = " << t0.next() << endl;
  }

  ballT D;
  switch (method) {
  case 0: {
    cout << "method = plain" << endl;
    auto support = vector<pointT>();
    if (serial)
      D = miniDiscPlainSerial(P, n, support, ballT());
    else
      D = miniDiscPlain(P, n, support, ballT());
    break;
  }
  case 1: {
    cout << "method = mtf" << endl;
    auto support = vector<pointT>();
    if (serial)
      D = miniDiscMtfSerial(P, n, support, ballT());
    else
      D = miniDiscMtf(P, n, support, ballT());
    break;
  }
  case 2: {
    cout << "method = mtf+pivot" << endl;
    auto support = vector<pointT>();
    if (serial)
      D = miniDiscPivotSerial(P, n, support, ballT());
    else
      D = miniDiscPivot(P, n, support, ballT());
    break;
  }
  default:
    cout << "invalid method, abort" << endl; abort();
  }

  cout << "seb-time = " << t0.stop() << endl;

  cout << D.radius() << ", center = " << D.center() << endl;

  cout << endl;

  check<dim,ballT>(&D, P, n);
}

template void miniDisc<2>(point<2>*, intT);
template void miniDisc<3>(point<3>*, intT);
template void miniDisc<4>(point<4>*, intT);
template void miniDisc<5>(point<5>*, intT);
template void miniDisc<6>(point<6>*, intT);
template void miniDisc<7>(point<7>*, intT);
template void miniDisc<8>(point<8>*, intT);
template void miniDisc<9>(point<9>*, intT);
