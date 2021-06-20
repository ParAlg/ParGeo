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
using namespace pargeo::graphIO;
using namespace pargeo::pointIO;

template<int dim>
void timeGraph(parlay::sequence<pargeo::point<dim>> &P, char const *outFile) {
  timer t; t.start();

  auto I = delaunayGraph<dim>(P);
  cout << "time = " << t.get_next() << endl;

  if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);

  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeGraph<2>(Points, oFile);
  } else {
    throw std::runtime_error("only supports 2d data set");
  }
}
