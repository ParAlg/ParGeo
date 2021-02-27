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
#include "parlay/hash_table.h"
#include "common/geometry.h"
#include <queue>
#include "hull.h"
#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

using namespace std;

floatT signedVolume(point<3> a, point<3> b, point<3> c, point<3> d) {
  return (a-d).dot(crossProduct3d(b-a, c-a))/6;
}

intT pivotOnEdge(intT q0, intT q1, parlay::sequence<point<3>> const &P) {
  intT p = 0;
  for (intT i=1; i<P.size(); ++i) {
    if (signedVolume(P[q0], P[q1], P[p], P[i]) > 0)
      p = i;
  }
  return p;
}

intT pivotOnEdge0(point<3> q0, point<3> q1, parlay::sequence<point<3>> const &P) {
  intT p = 0;
  for (intT i=1; i<P.size(); ++i) {
    if (signedVolume(q0, q1, P[p], P[i]) > 0)
      p = i;
  }
  return p;
}

parlay::sequence<facet3d> hull3d(parlay::sequence<point<3>> &P) {//todo change to hull 3d
  size_t n = P.size();

#ifdef WRITE
  ofstream myfile;
  myfile.open("point.txt", std::ofstream::trunc);
  for(auto p: P)
    myfile << p << endl;
  myfile.close();
#endif

  auto xCmp = [&](point<3> i, point<3> j) {
		     return i[0] < j[0];};
  intT p1 = parlay::min_element(P, xCmp) - &P[0];

  point<3> q = P[p1];
  q[1] += 1;//todo special case of having q in P where q only differ from p in x

#ifdef WRITE
  myfile.open("hull.txt", std::ofstream::trunc);
#endif

  auto H = parlay::sequence<facet3d>();

  // Find initial facet
  intT p2 = pivotOnEdge0(P[p1], q, P);
  intT p3 = pivotOnEdge(p1, p2, P);
  H.emplace_back(p1, p2, p3, P);

#ifdef WRITE
  myfile << P[p1] << endl;
  myfile << P[p2] << endl;
  myfile << P[p3] << endl;
#endif

  queue<pair<intT, intT>> Q;
  Q.push(make_pair(p2, p1));
  Q.push(make_pair(p3, p2));
  Q.push(make_pair(p1, p3));

  auto table = parlay::sequence<short>(n*n, 0);
  auto markProcessed = [&](intT p1, intT p2) {
		       table[n*p1+p2] = 1;
		     };
  auto processed = [&](intT p1, intT p2) {
		       return table[n*p1+p2];
		     };
  markProcessed(p1, p2);
  markProcessed(p2, p3);
  markProcessed(p3, p1);

  while (Q.size() > 0) {
    pair<intT, intT> e = Q.front();
    Q.pop();
    if (processed(e.first, e.second))
      continue;
    intT q = pivotOnEdge(e.first, e.second, P);
#ifdef WRITE
    myfile << P[q] << endl;
#endif
    H.emplace_back(e.first, q, e.second, P);
    Q.push(make_pair(H.back().b, H.back().a));
    Q.push(make_pair(H.back().c, H.back().b));
    Q.push(make_pair(H.back().a, H.back().c));
    markProcessed(H.back().a, H.back().b);
    markProcessed(H.back().b, H.back().c);
    markProcessed(H.back().c, H.back().a);
  }

#ifdef WRITE
  myfile.close();
#endif
  cout << "hull size = " << H.size() << endl;
  return H;
}

parlay::sequence<facet3d> hull3d(parlay::sequence<point<3>> &);
