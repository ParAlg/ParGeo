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
void timeGraph(parlay::sequence<pargeo::point<dim>> &P, double k, char const *outFile) {
  timer t; t.start();
  auto I = spanner<dim>(P, k);
  cout << "time = " << t.stop() << endl;
  if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-t <param>] [-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  double k = P.getOptionDoubleValue("-t",10); // spanner relaxation constant
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (k <= 1) throw std::runtime_error("Bad relaxation constant t");

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeGraph<2>(Points, k, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeGraph<3>(Points, k, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeGraph<4>(Points, k, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeGraph<5>(Points, k, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeGraph<6>(Points, k, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeGraph<7>(Points, k, oFile);
  } else {
    throw std::runtime_error("unsupported dimensionality");
  }
}
