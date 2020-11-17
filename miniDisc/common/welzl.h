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

#ifndef WELZL_H
#define WELZL_H

#include <vector>
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "prefix.h"
#include "miniDisc.h"
#include "geometry.h"

using namespace std;

template<int dim>
ball<dim> support2Ball(point<dim>* P, vector<point<dim>>& support) {
  typedef ball<dim> ballT;

  ballT B;
  if (B.isEmpty()) {
    if (support.size() == 0) {
      B = ballT(P, 2);
    } else if (support.size() == 1) {
      support.push_back(P[0]);
      B = ballT(&support[0], support.size());
      support.pop_back();
    } else { //>=2
      B = ballT(&support[0], support.size());
    }
  }
  return B;
}

template<int dim>
ball<dim> miniDiscPlainSerial(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
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
      B = miniDiscPlainSerial(P, i, support, B);
      support.pop_back();
    }
  }

  return B;
}

template<int dim>
ball<dim> miniDiscPlain(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B, intT* flag=NULL) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  if (n < 2000) return miniDiscPlainSerial(P, n, support, B);

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
                   B = miniDiscPlain(A, i, support, B, flag);
                   support.pop_back();
                 };
  parallel_prefix(P, n, process, cleanUp, freeFlag, flag);

  if(freeFlag) free(flag);
  return B;
}

#endif
