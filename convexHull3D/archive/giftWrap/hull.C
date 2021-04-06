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
#include "geometry/point.h"
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
template <class pt>
typename pt::floatT facetAngle(pt a1, pt b1, pt c1, pt c2) {
  pt v1 = crossProduct3d(b1-a1, c1-a1);
  pt v2 = crossProduct3d(b1-a1, c2-a1);
  typename pt::floatT dot = v1.dot(v2);
  typename pt::floatT angle = acos(dot / sqrt(v1.lenSqr()*v2.lenSqr()));
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
template <class pt>
typename pt::floatT pivotOnFacet0(pt p1, pt q1, pt q2, parlay::slice<pt*, pt*> P) {
  const pt* q = parlay::max_element(P, [&](pt a, pt b) {
					 auto ag1 = facetAngle(p1, q1, q2, a);
					 auto ag2 = facetAngle(p1, q1, q2, b);
					 return ag1 < ag2 ? true : false;
				       });
  pt qq = *q;
  return q-P.begin();
}

template <class pt>
size_t pivotOnFacet(size_t p1, size_t q1, size_t q2, parlay::slice<pt*, pt*> P) {
  const pt* q = parlay::max_element(P, [&](pt a, pt b) {
					 auto ag1 = facetAngle(P[p1], P[q1], P[q2], a);
					 auto ag2 = facetAngle(P[p1], P[q1], P[q2], b);
					 return ag1 < ag2 ? true : false;
				       });
  pt qq = *q;
  return q-P.begin();
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  using pt = pargeo::fpoint<3>;
  using facet3d = facet3d<pt>;

  size_t n = P.size();

#ifdef WRITE
  ofstream myfile;
  myfile.open("point.txt", ofstream::trunc);
  for(size_t p=0; p<P.size(); ++p)
    myfile << P[p] << p << endl;
  myfile.close();
#endif

  auto xCmp = [&](pt i, pt j) {
		return i[0] < j[0];};
  size_t p1 = parlay::min_element(P, xCmp) - &P[0];

  // todo special cases for colinear points
  pt q1 = P[p1];
  q1[1] += 1;
  pt q2 = P[p1];
  q2[2] += 1;
  size_t p2 = pivotOnFacet0(P[p1], q1, q2, make_slice(P));
  size_t p3 = pivotOnFacet0(P[p1], P[p2], q1, make_slice(P));

  auto H = parlay::sequence<facet3d>();
  H.emplace_back(p1, p2, p3, make_slice(P));

#ifdef WRITE
  myfile.open("hull.txt", ofstream::trunc);
  myfile << P[p1] << p1 << endl;
  myfile << P[p2] << p2 << endl;
  myfile << P[p3] << p3 << endl;
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
      size_t e1, e2, e3;
      if (j==0) {
	e2 = f->a; e1 = f->b; e3 = f->c;
      } else if (j==1) {
	e2 = f->b; e1 = f->c; e3 = f->a;
      } else {
	e2 = f->c; e1 = f->a; e3 = f->b;
      }

      if (T.processed(e1, e2))
	continue;
      size_t q = pivotOnFacet(e1, e2, e3, make_slice(P));
#ifdef WRITE
      myfile << P[q] << endl;
#endif
      H.emplace_back(e1, q, e2, make_slice(P));
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

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
