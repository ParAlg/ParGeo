#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "geometry/point.h"
#include "common/get_time.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "knnSearch/kdtKnn.h"

using namespace std;
using namespace benchIO;
using namespace pargeo;

template<int dim>
void timeKnn(parlay::sequence<pargeo::point<dim>> &P, size_t k, int rounds, char const *outFile) {
  timer t; t.start();
  parlay::sequence<size_t> I;
  for(int i=0; i<rounds; ++i) {
    I = kdtKnn::kdtKnn<dim, pargeo::point<dim>>(P, k+1);
    cout << "round-time = " << t.get_next() << endl;
  }
  t.stop();
  if (outFile != NULL) writeIntSeqToFile(I, outFile);
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
    timeKnn<2>(Points, k, rounds, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeKnn<3>(Points, k, rounds, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeKnn<4>(Points, k, rounds, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeKnn<5>(Points, k, rounds, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeKnn<6>(Points, k, rounds, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeKnn<7>(Points, k, rounds, oFile);
  } else {
    cout << "dim = " << dim << endl;
    throw std::runtime_error("dimension not supported");
  }
}
