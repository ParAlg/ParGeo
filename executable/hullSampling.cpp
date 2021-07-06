#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/sampling/sampling.h"

#include <iostream>
#include <algorithm>
#include "pargeo/parseCommandLine.h"
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "dataset/uniform.h"

using namespace pargeo;

int main(int argc, char* argv[]) {

  commandLine P(argc, argv, "[-n <data-size>] [-f <sample-fraction>]");
  std::cout << "[-n <data-size>] [-f <sample-fraction>]\n";

  size_t n = P.getOptionDoubleValue("-n", 1000000);
  double fraction = P.getOptionDoubleValue("-f", 0.001);

  auto pts = pargeo::uniformInPolyPoints<3, fpoint<3>>(n, 0);

  timer t; t.start();

  std::cout << "random\n";
  auto pRand = pargeo::hull3d::samplingHelper::randomSample<fpoint<3>>(parlay::make_slice(pts), fraction * n);
  std::cout << pRand.size() << "\n";
  std::cout << "h-fraction = " << pargeo::hull3d::sampling::test<fpoint<3>>(parlay::make_slice(pRand), 1) << "\n";
  std::cout << "time = " << t.get_next() << "\n";

  std::cout << "rand-proj\n";
  auto pRand2 = pargeo::hull3d::samplingHelper::randomProjection<fpoint<3>>(parlay::make_slice(pts), fraction * n);
  std::cout << pRand2.size() << "\n";
  std::cout << "h-fraction = " << pargeo::hull3d::sampling::test<fpoint<3>>(parlay::make_slice(pRand2), 1) << "\n";
  std::cout << "time = " << t.get_next() << "\n";

  std::cout << "grid-sample\n";
  auto pRand3 = pargeo::hull3d::samplingHelper::gridSample<fpoint<3>>(parlay::make_slice(pts), fraction * n);
  std::cout << pRand3.size() << "\n";
  std::cout << "h-fraction = " << pargeo::hull3d::sampling::test(parlay::make_slice(pRand3), 1) << "\n";
  std::cout << "time = " << t.get_next() << "\n";

  return 0;
}
