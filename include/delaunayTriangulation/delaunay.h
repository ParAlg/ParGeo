// The inteface for delaunay triangulation
#include "geometry.h"
#include "parlay/primitives.h"

namespace pbbsbench {

  using coord = double;
  using pointT = point2d<coord>;

  triangles<pointT> delaunay(parlay::sequence<pointT>& P);

}
