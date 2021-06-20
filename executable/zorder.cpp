#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/zorderSort.h"
#include "pargeo/parseCommandLine.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;

void timeSort2d(parlay::sequence<pargeo::point<2>> &P, char const *outFile) {
  timer t; t.start();
  auto sorted = pargeo::zorderSort2d<pargeo::point<2>>(P);
  cout << "time = " << t.stop() << endl;
  if (outFile != NULL) pargeo::pointIO::writePointsToFile(sorted, outFile);
}

void timeSort3d(parlay::sequence<pargeo::point<3>> &P, char const *outFile) {
  timer t; t.start();
  auto sorted = pargeo::zorderSort3d<pargeo::point<3>>(P);
  cout << "time = " << t.stop() << endl;
  if (outFile != NULL) pargeo::pointIO::writePointsToFile(sorted, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeSort2d(Points, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeSort3d(Points, oFile);
  } else {
    cout << "dim = " << dim << endl;
    throw std::runtime_error("dimension not supported");
  }
}
