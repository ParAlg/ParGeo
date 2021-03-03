// This code is part of the Pargeo Library
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

//#define WRITE // Write to file, visualize using python3 plot.py

#include <algorithm>
#include <queue>
#include <math.h>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "hull.h"
#include "pairHash.h"
#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

using namespace std;

/* Angle between facet(a1, c1, b1) & facet(a1, c2, b1)
 c1
   o
   | \
a1 o--o b1
    \ |
      o c2
 */
floatT facetAngle(point<3> a1, point<3> b1, point<3> c1, point<3> c2) {
  point<3> v1 = crossProduct3d(b1-a1, c1-a1);
  point<3> v2 = crossProduct3d(b1-a1, c2-a1);
  floatT dot = v1.dot(v2);
  floatT angle = acos(dot / sqrt(v1.lenSqr()*v2.lenSqr()));
  if (isnan(angle)) return 0.0;
  else return angle;
}

/* Pivot on p1-q1 of facet(p1, q2, q1)
 q2
  o
  | \
  o--o q1
 p1
 */
intT pivotOnFacet0(point<3> p1, point<3> q1, point<3> q2, parlay::sequence<point<3>> const &P) {
  typedef point<3> pt;

  const pt* q = parlay::max_element(P, [&](pt a, pt b) {
				   auto ag1 = facetAngle(p1, q1, q2, a);
				   auto ag2 = facetAngle(p1, q1, q2, b);
				   return ag1 < ag2 ? true : false;
				 });
  pt qq = *q;
  return q-P.begin();
}

intT pivotOnFacet(intT p1, intT q1, intT q2, parlay::sequence<point<3>> const &P) {
  typedef point<3> pt;

  const pt* q = parlay::max_element(P, [&](pt a, pt b) {
				   auto ag1 = facetAngle(P[p1], P[q1], P[q2], a);
				   auto ag2 = facetAngle(P[p1], P[q1], P[q2], b);
				   return ag1 < ag2 ? true : false;
				 });
  pt qq = *q;
  return q-P.begin();
}

parlay::sequence<facet3d> hull3d(parlay::sequence<point<3>> &P) {
  size_t n = P.size();

#ifdef WRITE
  ofstream myfile;
  myfile.open("point.txt", ofstream::trunc);
  for(auto p: P)
    myfile << p << endl;
  myfile.close();
#endif

  auto xCmp = [&](point<3> i, point<3> j) {
		     return i[0] < j[0];};
  intT p1 = parlay::min_element(P, xCmp) - &P[0];

  // todo special cases for colinear points
  point<3> q1 = P[p1];
  q1[1] += 1;
  point<3> q2 = P[p1];
  q2[2] += 1;
  intT p2 = pivotOnFacet0(P[p1], q1, q2, P);
  intT p3 = pivotOnFacet0(P[p1], P[p2], q1, P);

  auto H = parlay::sequence<facet3d>();
  H.emplace_back(p1, p2, p3, P);

#ifdef WRITE
  myfile.open("hull.txt", ofstream::trunc);
  myfile << P[p1] << endl;
  myfile << P[p2] << endl;
  myfile << P[p3] << endl;
#endif

  queue<facet3d*> Q;
  Q.push(&H[0]);

  auto T = pairHash(4*n);

  // p1, p2, p3 are not oriented
  T.mark(p2, p1);
  T.mark(p3, p2);
  T.mark(p1, p3);

  while (Q.size() > 0) {
    facet3d *f = Q.front();
    Q.pop();
    for (int j=0; j<3; ++j) {
      intT e1, e2, e3;
      if (j==0) {
	e2 = f->a; e1 = f->b; e3 = f->c;
      } else if (j==1) {
	e2 = f->b; e1 = f->c; e3 = f->a;
      } else {
	e2 = f->c; e1 = f->a; e3 = f->b;
      }

      if (T.processed(e1, e2))
	continue;
      intT q = pivotOnFacet(e1, e2, e3, P);
#ifdef WRITE
      myfile << P[q] << endl;
#endif
      H.emplace_back(e1, q, e2, P);
      Q.push(&H[H.size()-1]);
      T.mark(H.back().a, H.back().b);
      T.mark(H.back().b, H.back().c);
      T.mark(H.back().c, H.back().a);
    }
  }

#ifdef WRITE
  myfile.close();
#endif
  cout << "hull size = " << H.size() << endl;
  return H;
}

parlay::sequence<facet3d> hull3d(parlay::sequence<point<3>> &);
