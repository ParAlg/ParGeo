// The inteface for delaunay triangulation
#include "common/geometry.h"
#include "parlay/primitives.h"

using coord = double;
using pointT = point2d<coord>;

triangles<pointT> delaunay(parlay::sequence<pointT>& P);
