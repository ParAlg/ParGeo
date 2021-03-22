#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "incremental.h"
#include "hull.h"

using namespace std;

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &P) {
  size_t n = P.size();

  //auto H = incrementHull3dSerial<pargeo::fpoint<3>>(make_slice(P));
  auto H = incrementHull3d<pargeo::fpoint<3>>(make_slice(P));

  //cout << "hull size = " << H.size() << endl;
  return H;
}

parlay::sequence<facet3d<pargeo::fpoint<3>>> hull3d(parlay::sequence<pargeo::fpoint<3>> &);
