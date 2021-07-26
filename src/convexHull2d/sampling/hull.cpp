#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"

#include "convexHull2d/quickHull/hull.h"
#include "convexHull2d/sampling/hull.h"
#include "convexHull2d/sampling/sampling.h"
#include "convexHull2d/sampling/filter.h"

// #define SAMPLE_HULL_VERBOSE

template<class pointT>
parlay::sequence<size_t>
pargeo::hull2d::sampling::compute(parlay::slice<pointT*, pointT*> P, double fraction) {

  // auto tmp = pargeo::hull2d::quickHull::compute(P);
  // for (auto p: tmp) {
  //   std::cout << P[p] << "\n";
  // }
  // std::cout << "\n";
  // abort();

  if (P.size() < 1000)
    return pargeo::hull2d::quickHull::compute(P);

#ifdef SAMPLE_HULL_VERBOSE
  pargeo::timer t; t.start();
#endif

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "sample-size = " << sampleSize << "\n";
#endif

  auto sample = pargeo::hull2d::samplingHelper::randomProjection<pointT>(P, sampleSize);

  parlay::sequence<size_t> sampleHull =
    pargeo::hull2d::quickHull::compute(parlay::make_slice(sample));

  // for (auto p: sampleHull) {
  //   std::cout << P[p] << "\n";
  // }
  // std::cout << "\n";
  // abort();

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "precompute-time = " << t.get_next() << "\n";
  std::cout << "h = " << sampleHull.size() << "\n";
#endif

  auto remain =
    pargeo::hull2d::samplingHelper::filter2<pointT>(P, parlay::make_slice(sample), parlay::make_slice(sampleHull));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "filter-time = " << t.get_next() << "\n";
  std::cout << "remain = " << remain.size() << "\n";
  std::cout << "fraction = " << double(remain.size())/P.size() << "\n";
#endif

  auto hull = pargeo::hull2d::quickHull::compute(parlay::make_slice(remain));

#ifdef SAMPLE_HULL_VERBOSE
  std::cout << "final-hull-time = " << t.stop() << "\n";
#endif
  return hull; // todo note the returned result index is wrt the remaining points ("remain")
}

template parlay::sequence<size_t>
pargeo::hull2d::sampling::compute(parlay::slice<pargeo::fpoint<2>*, pargeo::fpoint<2>*>, double);

template parlay::sequence<size_t>
pargeo::hull2d::sampling::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*>, double);
