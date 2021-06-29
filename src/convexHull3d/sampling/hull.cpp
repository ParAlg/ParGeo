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

#include "convexHull3d/vertex.h"
#include "convexHull3d/serialHull.h"
#include "convexHull3d/samplingHull.h"
#include "convexHull3d/internal/sampling.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace pargeo;

// #define SAMPLE_HULL_VERBOSE

parlay::sequence<facet3d<pargeo::fpoint<3>>>
pargeo::hull3dSampling(parlay::sequence<pargeo::fpoint<3>> &P, double fraction) {
  using namespace std;
  using namespace parlay;
  using pt = pargeo::hullInternal::vertex;
  using floatT = typename pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;

  if (P.size() < 1000) return hull3dSerial(P);

#ifdef SAMPLE_HULL_VERBOSE
  timer t; t.start();
#endif

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

  auto sample = pargeo::hullInternal::randomSample(P, sampleSize);

  parlay::sequence<facet3d<pointT>> sampleHull = hull3dSerial(sample);

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "precompute-time = " << t.get_next() << "\n";
  std::cout << "h = " << sampleHull.size() << "\n";
#endif

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

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "filter-time = " << t.get_next() << "\n";
  std::cout << "remain = " << remain.size() << "\n";
  std::cout << "fraction = " << double(remain.size())/P.size() << "\n";
#endif

  auto hull = hull3dSerial(remain);

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "final-hull-time = " << t.stop() << "\n";
#endif
  return hull;
}
