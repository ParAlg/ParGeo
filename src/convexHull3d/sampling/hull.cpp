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
#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/sampling/sampling.h"
#include "convexHull3d/sampling/filter.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"

#include "pargeo/getTime.h"
#include "pargeo/point.h"

// #define SAMPLE_HULL_VERBOSE

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::sampling::compute(parlay::slice<pointT*, pointT*> P, double fraction) {

  if (P.size() < 1000)
    return pargeo::hull3d::serialQuickHull::compute(P);

#ifdef SAMPLE_HULL_VERBOSE
  pargeo::timer t; t.start();
#endif

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

  auto sample = pargeo::hull3d::samplingHelper::randomProjection<pointT>(P, sampleSize);

  parlay::sequence<pargeo::hull3d::facet<pointT>> sampleHull =
    pargeo::hull3d::serialQuickHull::compute(make_slice(sample));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "precompute-time = " << t.get_next() << "\n";
  std::cout << "h = " << sampleHull.size() << "\n";
#endif

  auto remain =
    pargeo::hull3d::samplingHelper::filter2<pointT>(P, parlay::make_slice(sampleHull));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "filter-time = " << t.get_next() << "\n";
  std::cout << "remain = " << remain.size() << "\n";
  std::cout << "fraction = " << double(remain.size())/P.size() << "\n";
#endif

  auto hull = pargeo::hull3d::serialQuickHull::compute(make_slice(remain));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "final-hull-time = " << t.stop() << "\n";
#endif
  return hull;
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::sampling::compute(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*> P, double fraction);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::sampling::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*> P, double fraction);



template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::sampling::random(parlay::slice<pointT*, pointT*> P, double fraction) {
  if (P.size() < 1000) return pargeo::hull3d::serialQuickHull::compute<pointT>(P);
  size_t sampleSize = std::max(size_t(P.size() * fraction), size_t(5));
  auto sample = pargeo::hull3d::samplingHelper::randomSample<pointT>(P, sampleSize);
  parlay::sequence<pargeo::hull3d::facet<pointT>> sampleHull =
    pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(sample));
  auto remain = pargeo::hull3d::samplingHelper::filter2(P, make_slice(sampleHull));
  return pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(remain));
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::sampling::random(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, double);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::sampling::random(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, double);

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::sampling::grid(parlay::slice<pointT*, pointT*> P, double fraction) {
  if (P.size() < 1000) return pargeo::hull3d::serialQuickHull::compute<pointT>(P);
  size_t sampleSize = std::max(size_t(P.size() * fraction), size_t(5));
  auto sample = pargeo::hull3d::samplingHelper::gridSample<pointT>(P, sampleSize);
  parlay::sequence<pargeo::hull3d::facet<pointT>> sampleHull =
    pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(sample));
  auto remain = pargeo::hull3d::samplingHelper::filter2(P, make_slice(sampleHull));
  return pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(remain));
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::sampling::grid(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, double);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::sampling::grid(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, double);

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::sampling::projection(parlay::slice<pointT*, pointT*> P, double fraction) {
  if (P.size() < 1000) return pargeo::hull3d::serialQuickHull::compute<pointT>(P);
  size_t sampleSize = std::max(size_t(P.size() * fraction), size_t(5));
  auto sample = pargeo::hull3d::samplingHelper::randomProjection<pointT>(P, sampleSize);
  parlay::sequence<pargeo::hull3d::facet<pointT>> sampleHull =
    pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(sample));
  auto remain = pargeo::hull3d::samplingHelper::filter2(P, make_slice(sampleHull));
  return pargeo::hull3d::serialQuickHull::compute<pointT>(make_slice(remain));
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::sampling::projection(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, double);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::sampling::projection(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, double);
