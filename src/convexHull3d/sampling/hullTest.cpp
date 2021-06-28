#include "convexHull3d/serialHull.h"
#include "convexHull3d/samplingHull.h"
#include "convexHull3d/vertex.h"

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

using namespace pargeo;

double pargeo::testHull(parlay::sequence<pargeo::fpoint<3>> &P, double fraction) {
  using namespace parlay;
  using pt = pargeo::hullInternal::vertex;

  if (P.size() < 5) return 1.0;

  size_t sampleSize = P.size() * fraction;
  sampleSize = std::max(sampleSize, size_t(5));

  auto sample = parlay::tabulate(sampleSize, [&](size_t i) {
					       return P[parlay::hash64(i) % P.size()];
					     });

  sequence<pargeo::fpoint<3>> hullVertices = hullInternal::hull3dSerialInternal3(sample);

  return double(hullVertices.size()) / double(sampleSize);
}
