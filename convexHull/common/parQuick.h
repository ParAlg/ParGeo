// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#ifndef PAR_QUICK_H
#define PAR_QUICK_H

#include <algorithm>
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "serialQuick.h"
using namespace std;
using namespace sequence;

struct triangArea {
  intT l, r;
  point2d* P;
  intT* I;
  triangArea(intT* _I, point2d* _P, intT _l, intT _r) : I(_I), P(_P), l(_l), r(_r) {}
  double operator() (intT i) {return triArea(P[l], P[r], P[I[i]]);}
};

intT quickHull(intT* I, intT* Itmp, point2d* P, intT n, intT l, intT r, intT depth) {
  if (n < 10000 || depth == 0)
    return serialQuickHullHelper(I, P, n, l, r);
  else {

    intT idx = maxIndex<double>((intT)0,n,greater<double>(),triangArea(I,P,l,r));
    intT maxP = I[idx];

    intT n1 = filter(I, Itmp,    n, aboveLine(P, l, maxP));
    intT n2 = filter(I, Itmp+n1, n, aboveLine(P, maxP, r));

    intT m1, m2;
    par_do([&](){m1 = quickHull(Itmp, I ,P, n1, l, maxP, depth-1);},
	   [&](){m2 = quickHull(Itmp+n1, I+n1, P, n2, maxP, r, depth-1);});

    parallel_for (0, m1, [&](intT i) {I[i] = Itmp[i];});
    I[m1] = maxP;
    parallel_for (0, m2, [&](intT i) {I[i+m1+1] = Itmp[i+n1];});
    return m1+1+m2;
  }
}

struct makePair {
  pair<intT,intT> operator () (intT i) { return pair<intT,intT>(i,i);}
};

struct minMaxIndex {
  point2d* P;
  minMaxIndex (point2d* _P) : P(_P) {}
  pair<intT,intT> operator () (pair<intT,intT> l, pair<intT,intT> r) {
    intT minIndex =
      (P[l.first].x() < P[r.first].x()) ? l.first :
      (P[l.first].x() > P[r.first].x()) ? r.first :
      (P[l.first].y() < P[r.first].y()) ? l.first : r.first;
    intT maxIndex = (P[l.second].x() > P[r.second].x()) ? l.second : r.second;
    return pair<intT,intT>(minIndex, maxIndex);
  }
};

_seq<intT> hullInternal(point2d* P, intT n) {
  static const intT DEPTH = 10;

  pair<intT,intT> minMax = reduce<pair<intT,intT> >((intT)0,n,minMaxIndex(P), makePair());
  intT l = minMax.first;
  intT r = minMax.second;
  bool* fTop = newA(bool,n);
  bool* fBot = newA(bool,n);
  intT* I = newA(intT, n);
  intT* Itmp = newA(intT, n);
  parallel_for(0, n, [&](intT i) {
		       Itmp[i] = i;
		       double a = triArea(P[l],P[r],P[i]);
		       fTop[i] = a > 0;
		       fBot[i] = a < 0;
		     });

  intT n1 = pack(Itmp, I, fTop, n);
  intT n2 = pack(Itmp, I+n1, fBot, n);
  free(fTop); free(fBot);

  intT m1; intT m2;
  par_do([&](){m1 = quickHull(I, Itmp, P, n1, l, r, DEPTH);},
	 [&](){m2 = quickHull(I+n1, Itmp+n1, P, n2, r, l, DEPTH);});

  parallel_for (0, m1, [&](intT i) {Itmp[i+1] = I[i];});
  parallel_for (0, m2, [&](intT i) { Itmp[i+m1+2] = I[i+n1];});
  free(I);
  Itmp[0] = l;
  Itmp[m1+1] = r;

  return _seq<intT>(Itmp, m1+2+m2);
}

#endif
