#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/bruteforce/hull.h"
#include "convexHull3d/pseudo/hull.h"

#include <iostream>
#include <algorithm>
#include "pargeo/parseCommandLine.h"
#include "parlay/parallel.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"

using namespace pargeo;
using namespace pargeo::pointIO;

template <class pt>
void timeHull(parlay::sequence<pt> &P, int rounds, char const *outFile) {
  timer t; t.start();

  auto H = pargeo::hull3d::pseudo::compute<pt>(make_slice(P));

  std::cout << "hull-time = " << t.get_next() << "\n";
  std::cout << "hull-size = " << H.size() << "\n";

  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",0);

  int dim = readHeader(iFile);
  if (dim != 3) {
    std::cout << "Error, convexHull3D only takes 3d inputs, abort.\n";
    abort();
  }

  // todo precision
  parlay::sequence<pargeo::fpoint<3>> Points = readPointsFromFile<pargeo::fpoint<3>>(iFile);
  timeHull<pargeo::fpoint<3>>(Points, rounds, oFile);
  // parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
  // timeHull<pargeo::point<3>>(Points, rounds, oFile);
}
