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

void timeSort2d(parlay::sequence<pargeo::point<2>> &P, size_t k, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    auto I = pargeo::zorderSort2d<pargeo::point<2>>(P);
    // for (auto x: I)
    //   cout << x << endl;
    cout << "zorder-time = " << t.get_next() << endl;
  }
  t.stop();
}

void timeSort3d(parlay::sequence<pargeo::point<3>> &P, size_t k, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    auto I = pargeo::zorderSort3d<pargeo::point<3>>(P);
    // for (auto x: I)
    //   cout << x << endl;
    cout << "zorder-time = " << t.get_next() << endl;
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
    timeSort2d(Points, k, rounds, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeSort3d(Points, k, rounds, oFile);
  } else {
    cout << "dim = " << dim << endl;
    throw std::runtime_error("dimension not supported");
  }
}
