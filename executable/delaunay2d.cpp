#include <iostream>
#include <algorithm>
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"
#include "delaunayTriangulation/geometry.h"
#include "delaunayTriangulation/delaunay.h"
#include "parlay/primitives.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;
using namespace pbbsbench;

int main(int argc, char* argv[]) {
  using pointT = pbbsbench::point2d<double>;

  commandLine C(argc,argv,"./delaunay <inFile>");
  char* iFile = C.getArgument(0);

  parlay::sequence<pointT> P = readPointsFromFile<pointT>(iFile);
  triangles<pointT> R;

  timer t; t.start();
  R = delaunay(P);
  cout << "time = " << t.stop() << endl;
}
