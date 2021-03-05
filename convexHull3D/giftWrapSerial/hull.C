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
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pairHash.h"
#include "geometry/point.h"
#include "geometry/algebra.h"
#include <queue>
#include "hull.h"
#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

using namespace std;

/* Signed volume of an oriented tetrahedron (example below is positive).
    d
    o
   /|\
  / | o b
 o--o/
 a   c
 */
template <class pt>
typename pt::floatT signedVolume(pt a, pt b, pt c, pt d) {
  return (a-d).dot(crossProduct3d(b-a, c-a))/6;
}

template <class pt>
size_t pivotOnEdge(size_t q0, size_t q1, parlay::slice<pt*, pt*> P) {
  size_t p = 0;
  for (size_t i=1; i<P.size(); ++i) {
    if (signedVolume(P[q0], P[q1], P[p], P[i]) > 1e-7) //numeric knob
      p = i;
  }
  return p;
}

template <class pt>
size_t pivotOnEdge0(pt q0, pt q1, parlay::slice<pt*, pt*> P) {
  size_t p = 0;
  for (size_t i=1; i<P.size(); ++i) {
    if (signedVolume(q0, q1, P[p], P[i]) > 1e-7) //numeric knob
      p = i;
  }
  return p;
}

parlay::sequence<facet3d<pargeo::point<3>>> hull3d(parlay::sequence<pargeo::point<3>> &P) {
  using pt = pargeo::point<3>;
  using facet3d = facet3d<pt>;

  size_t n = P.size();

#ifdef WRITE
  ofstream myfile;
  myfile.open("point.txt", std::ofstream::trunc);
  for(size_t p=0; p<P.size(); ++p)
    myfile << P[p] << p << endl;
  myfile.close();
#endif

  auto xCmp = [&](pt i, pt j) {
		return i[0] < j[0];};
  size_t p1 = parlay::min_element(P, xCmp) - &P[0];

  pt q = P[p1];
  q[1] += 1;//todo special case of having q in P where q only differ from p in x

#ifdef WRITE
  myfile.open("hull.txt", std::ofstream::trunc);
#endif

  auto H = parlay::sequence<facet3d>();

  // Find initial facet
  size_t p2 = pivotOnEdge0(P[p1], q, make_slice(P));
  size_t p3 = pivotOnEdge(p1, p2, make_slice(P));
  H.emplace_back(p1, p2, p3, make_slice(P));

#ifdef WRITE
  myfile << P[p1] << p1 << endl;
  myfile << P[p2] << p1 << endl;
  myfile << P[p3] << p1 << endl;
#endif

  queue<pair<size_t, size_t>> Q;
  Q.push(make_pair(p2, p1));
  Q.push(make_pair(p3, p2));
  Q.push(make_pair(p1, p3));

  auto T = pairHash(4*n);
  T.mark(p1, p2);
  T.mark(p2, p3);
  T.mark(p3, p1);

  while (Q.size() > 0) {
    pair<size_t, size_t> e = Q.front();
    Q.pop();
    if (T.processed(e.first, e.second))
      continue;
    size_t q = pivotOnEdge(e.first, e.second, make_slice(P));
#ifdef WRITE
    myfile << P[q] << q << endl;
#endif
    H.emplace_back(e.first, q, e.second, make_slice(P));
    Q.push(make_pair(H.back().b, H.back().a));
    Q.push(make_pair(H.back().c, H.back().b));
    Q.push(make_pair(H.back().a, H.back().c));
    T.mark(H.back().a, H.back().b);
    T.mark(H.back().b, H.back().c);
    T.mark(H.back().c, H.back().a);
  }

#ifdef WRITE
  myfile.close();
#endif
  cout << "hull size = " << H.size() << endl;
  return H;
}

parlay::sequence<facet3d<pargeo::point<3>>> hull3d(parlay::sequence<pargeo::point<3>> &);
