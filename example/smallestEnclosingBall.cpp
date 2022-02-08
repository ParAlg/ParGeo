#include "enclosingBall/sampling/seb.h"

#include <algorithm>
#include <iostream>
#include "parlay/parallel.h"
#include "pargeo/parseCommandLine.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"

template <int dim>
void callSeb(parlay::sequence<pargeo::point<dim>> &P) {

  auto D = pargeo::seb::sampling::compute<dim>(make_slice(P));

  std::cout << D.radius() << ", center = " << D.center() << "\n";
}

int main(int argc, char* argv[]) {
  pargeo::commandLine P(argc, argv, "<inFile>");
  char* iFile = P.getArgument(0);

  int dim = pargeo::pointIO::readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<2>>(iFile);
    callSeb<2>(P);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(iFile);
    callSeb<3>(P);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<4>>(iFile);
    callSeb<4>(P);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<5>>(iFile);
    callSeb<5>(P);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<6>>(iFile);
    callSeb<6>(P);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<7>>(iFile);
    callSeb<7>(P);
  } else {
    std::cout << "dim = " << dim << "\n";
    throw std::runtime_error("Dimension not supported.");
  }
}
