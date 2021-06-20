#include "convexHull3d/hull.h"

#include <iostream>
#include <algorithm>
#include "pargeo/parseCommandLine.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/point.h"
#include "parlay/parallel.h"

using namespace pargeo;
using namespace pargeo::pointIO;

template <class pt>
void timeHull(parlay::sequence<pt> &P, double fraction, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    auto H = hull3dSampling(P, fraction);
    std::cout << "hull-size = " << H.size() << "\n";
    std::cout << "round-time = " << t.get_next() << "\n";
  }
  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  double fraction = P.getOptionDoubleValue("-f",0);
  int rounds = P.getOptionIntValue("-r",1);

  int dim = readHeader(iFile);
  if (dim != 3) {
    throw std::runtime_error("Error, convexHull3D only takes 3d inputs");
  }

  parlay::sequence<pargeo::fpoint<3>> Points = readPointsFromFile<pargeo::fpoint<3>>(iFile);
  timeHull<pargeo::fpoint<3>>(Points, fraction, rounds, oFile);
}
