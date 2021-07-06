#include "convexHull3d/gift/hull.h"
#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/bruteforce/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/pseudo/hull.h"
#include "convexHull3d/parallelQuickHull/hull.h"
#include "convexHull3d/divideConquer/hull.h"

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

  //auto H = pargeo::hull3d::parallelQuickHull::compute<pt>(make_slice(P));
  //auto H = pargeo::hull3d::divideConquer::compute<pt>(make_slice(P));
  //auto H = pargeo::hull3d::sampling::compute<pt>(make_slice(P));
  auto H = pargeo::hull3d::pseudo::compute<pt>(make_slice(P));
  //auto H = pargeo::hull3d::gift::compute<pt>(make_slice(P));
  //auto H = pargeo::hull3d::serialQuickHull::compute<pt>(make_slice(P));
  //auto H = pargeo::hull3d::bruteforce::compute(make_slice(P));

  std::cout << "hull-time = " << t.get_next() << "\n";
  std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dGift(P);
  // std::cout << "serial-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dIncremental(P);
  // std::cout << "incremental-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dConcurrent(P);
  // std::cout << "concurrent-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dPseudo(P);
  // std::cout << "pseudo-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dGrid(P);
  // std::cout << "grid-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // H = hull3dSearch(P);
  // std::cout << "search-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // auto H = hull3dSampling(P, 0.01);
  // std::cout << "sampling-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";
  // auto H = hull3dPseudo(P);
  // auto H = hull3dParallel(P, rounds);
  //auto H = hull3dIncremental(P, rounds);
  // auto H = hull3dPseudo(P);
  // std::cout << "hull-time = " << t.get_next() << "\n";
  // std::cout << "hull-size = " << H.size() << "\n";

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

  parlay::sequence<pargeo::fpoint<3>> Points = readPointsFromFile<pargeo::fpoint<3>>(iFile);
  timeHull<pargeo::fpoint<3>>(Points, rounds, oFile);
  // parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
  // timeHull<pargeo::point<3>>(Points, rounds, oFile);
}
