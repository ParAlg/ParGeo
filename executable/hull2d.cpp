#include "convexHull2d/quickHull/hull.h"
#include "convexHull2d/bruteforce/hull.h"
#include "convexHull2d/divideConquer/hull.h"

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

  //auto H = pargeo::hull2d::quickHull::compute<pt>(make_slice(P));
  //auto H = pargeo::hull2d::bruteforce::compute<pt>(make_slice(P));
  auto H = pargeo::hull2d::divideConquer::compute<pt>(make_slice(P));

  // for (auto h: H) {
  //   std::cout << h << " ";
  // }
  // std::cout << "\n";

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
  if (dim != 2) {
    std::cout << "Error, convexHull2D only takes 2d inputs, abort.\n";
    abort();
  }

  parlay::sequence<pargeo::fpoint<2>> Points = readPointsFromFile<pargeo::fpoint<2>>(iFile);
  timeHull<pargeo::fpoint<2>>(Points, rounds, oFile);
}
