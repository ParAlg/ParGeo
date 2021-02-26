#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "common/geometry.h"
#include "hull.h"

using namespace std;

void hull(parlay::sequence<point<3>> const &P) {//todo change to hull 3d
  typedef int intT;
  typedef double floatT;
  size_t n = P.size();

  auto xCmp = [&](point<3> i, point<3> j) {
		     return i[0] < j[0];};

  intT xMin = parlay::min_element(P, xCmp) - &P[0];
  cout << "xMin = " << P[xMin] << endl;
}

void hull(parlay::sequence<point<3>> const &);
