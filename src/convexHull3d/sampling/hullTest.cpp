#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/vertex.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

template<class pointT>
double
pargeo::hull3d::sampling::test(parlay::slice<pointT*, pointT*> P, double fraction) {
  using namespace parlay;

  if (P.size() < 5) return 1.0;

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

  auto sample = parlay::tabulate(sampleSize, [&](size_t i) {
					       return P[parlay::hash64(i) % P.size()];
					     });

  auto hullVertices = pargeo::hull3d::serialQuickHull::computeVertex<pointT>(make_slice(sample));

  return double(hullVertices.size()) / double(sampleSize);
}

template double
pargeo::hull3d::sampling::test(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>, double);

template double
pargeo::hull3d::sampling::test(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>, double);
