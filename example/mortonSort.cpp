#include <iostream>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

#include "mortonSort/mortonSort.h"

void sort2d(parlay::sequence<pargeo::point<2>> &P, char const *outFile) {
  auto sorted = pargeo::zorderSort2d<pargeo::point<2>>(P);
  if (outFile != NULL) pargeo::pointIO::writePointsToFile(sorted, outFile);
}

void sort3d(parlay::sequence<pargeo::point<3>> &P, char const *outFile) {
  auto sorted = pargeo::zorderSort3d<pargeo::point<3>>(P);
  if (outFile != NULL) pargeo::pointIO::writePointsToFile(sorted, outFile);
}

int main(int argc, char* argv[]) {
  pargeo::commandLine P(argc,argv,"[-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  int dim = pargeo::pointIO::readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<2>>(iFile);
    sort2d(P, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(iFile);
    sort3d(P, oFile);
  } else {
    std::cout << "dim = " << dim << std::endl;
    throw std::runtime_error("Dimension not supported.");
  }
}
