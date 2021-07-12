#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/scan/seb.h"

#include <iostream>
#include <algorithm>
#include "pargeo/parseCommandLine.h"
#include "parlay/parallel.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"

using namespace pargeo;
using namespace pargeo::pointIO;

template <int dim>
void timeSeb(parlay::sequence<pargeo::point<dim>> &P, char const *outFile) {
  timer t; t.start();

  //auto D = pargeo::seb::welzlMtfPivot::compute<dim>(make_slice(P));
  auto D = pargeo::seb::scan::compute<dim>(make_slice(P));
  std::cout << D.radius() << ", center = " << D.center() << "\n";
  std::cout << "seb-time = " << t.get_next() << "\n";
  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeSeb<2>(Points, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeSeb<3>(Points, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeSeb<4>(Points, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeSeb<5>(Points, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeSeb<6>(Points, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeSeb<7>(Points, oFile);
  } else {
    std::cout << "dim = " << dim << "\n";
    throw std::runtime_error("dimension not supported");
  }
}
