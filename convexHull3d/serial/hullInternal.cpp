#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "convexHull3d/incremental.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/pointVertex.h"
#include "convexHull3d/hull.h"

parlay::sequence<pointVertex> hull3dInternalSerial(parlay::sequence<pointVertex> &Q) {
  using namespace std;
  using namespace parlay;
  using facetT = facet3d<pargeo::fpoint<3>>;
  using vertexT = pointVertex;

  // Create an initial simplex
  auto origin = pointOrigin();
  auto linkedHull = new _hull<linkedFacet3d<vertexT>, vertexT, pointOrigin>(make_slice(Q), origin, true);

  incrementHull3dSerial<linkedFacet3d<vertexT>, vertexT>(linkedHull);

  // linkedHull is translated
  // getHull will undo the translation
  auto vertices = linkedHull->getHullVertices<vertexT>();
  delete linkedHull;
  return vertices;
}

parlay::sequence<pointVertex> hull3dInternalSerial(parlay::sequence<pointVertex> &);
