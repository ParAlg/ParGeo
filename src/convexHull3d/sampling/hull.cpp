// This code is part of the Pargeo Library
// Copyright (c) 2021 Yiqiu Wang and the Pargeo Team
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

#include "convexHull3d/hull.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace pargeo;

parlay::sequence<facet3d<pargeo::fpoint<3>>>
pargeo::hull3dSampling(parlay::sequence<pargeo::fpoint<3>> &P, double fraction) {
  using namespace std;
  using namespace parlay;
  using pt = vertex;
  using floatT = typename pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;

  if (P.size() < 1000) return hull3dSerial(P);

  timer t; t.start();

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

  auto sample = parlay::tabulate(sampleSize, [&](size_t i) {
					       return P[parlay::hash64(i) % P.size()];
					     });

  parlay::sequence<facet3d<pointT>> sampleHull = hull3dSerial(sample);
  std::cout << "precompute-time = " << t.get_next() << "\n";
  std::cout << "h = " << sampleHull.size() << "\n";

  auto interiorPt = (sampleHull[0].a + sampleHull[sampleSize % sampleHull.size()].a) / 2;

  sequence<pointT> area(sampleHull.size());
  parallel_for(0, sampleHull.size(),
	       [&](size_t i){
		 auto f = sampleHull[i];
		 area[i] = pargeo::crossProduct3d(f.b-f.a, f.c-f.a);
	       });

  auto isOut = [&](pointT p) {
		 for (size_t i = 0; i < sampleHull.size(); ++i) {
		   auto f = sampleHull[i];
		   if (((f.a - interiorPt) - (p - interiorPt)).dot(area[i]) > 0) {
		     return true;
		   }
		 }
		 return false;
	       };

  auto remain = parlay::filter(make_slice(P), isOut);
  std::cout << "filter-time = " << t.get_next() << "\n";

  std::cout << "f = " << double(remain.size())/P.size() << "\n";

  auto hull = hull3dSerial(remain);
  std::cout << "final-hull-time = " << t.stop() << "\n";

  return hull;
}
