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
void timeGraph(parlay::sequence<pargeo::point<dim>> &P, double b, char const *outFile) {
  timer t; t.start();
  auto I = pargeo::betaSkeleton<dim>(P, b);
  cout << "time = " << t.stop() << endl;
  if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-b <param>] [-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  double b = P.getOptionDoubleValue("-b",1.0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeGraph<2>(Points, b, oFile);
  } else {
    throw std::runtime_error("unsupported dimensionality");
  }
}
