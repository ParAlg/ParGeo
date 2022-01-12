#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

template<int dim>
void knn(parlay::sequence<pargeo::point<dim>> &P, size_t k, char const *outFile) {
  //timer t; t.start();
  parlay::sequence<size_t> I =
    pargeo::kdTreeKnn<dim, pargeo::point<dim>>(P, k);
  //cout << "time = " << t.stop() << endl;
  if (outFile != NULL) pargeo::IO::writeIntSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  pargeo::commandLine P(argc,argv,"[-k <param>] [-o <outFile>] <inFile>\n k = 1 will just return self.");
  char* iFile = P.getArgument(0);
  size_t k = P.getOptionIntValue("-k",1);
  char* oFile = P.getOptionValue("-o");

  int dim = pargeo::pointIO::readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<2>>(iFile);
    knn<2>(P, k, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(iFile);
    knn<3>(P, k, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<4>>(iFile);
    knn<4>(P, k, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<5>>(iFile);
    knn<5>(P, k, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<6>>(iFile);
    knn<6>(P, k, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<7>>(iFile);
    knn<7>(P, k, oFile);
  } else {
    std::cout << "dim = " << dim << std::endl;
    throw std::runtime_error("Dimension not supported.");
  }
}
