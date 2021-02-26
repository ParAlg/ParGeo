#include "common/geometry.h"
#include "parlay/primitives.h"

typedef int intT;
typedef double floatT;

struct facet3d {
  // Indices into global point<3>*
  // Clockwise order
  intT a, b, c;
};

void hull(parlay::sequence<point<3>> const &S);
