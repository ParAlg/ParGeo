#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/graphIO.h"
#include "pargeo/parseCommandLine.h"
#include "spatialGraph/spatialGraph.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;
using namespace pargeo::graphIO;

template<int dim>
void timeGraph(parlay::sequence<pargeo::point<dim>> &P, size_t k, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    auto I = knnGraph<dim>(P, k);
    cout << "round-time = " << t.get_next() << endl;
    if (i == rounds-1 && outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
  }
  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  size_t k = P.getOptionIntValue("-k",1);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeGraph<2>(Points, k, rounds, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeGraph<3>(Points, k, rounds, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeGraph<4>(Points, k, rounds, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeGraph<5>(Points, k, rounds, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeGraph<6>(Points, k, rounds, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeGraph<7>(Points, k, rounds, oFile);
  } else {
    throw std::runtime_error("unsupported dimensionality");
  }
}
