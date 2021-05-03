#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/kdTreeRange.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;

template<int dim>
void timeKnn(parlay::sequence<pargeo::point<dim>> &P, double k, int rounds, char const *outFile) {
  timer t; t.start();

  for(int i=0; i<rounds; ++i) {

    kdNode<dim, point<dim>>* tree = buildKdt<dim, point<dim>>(P, true);

    cout << "build-time = " << t.get_next() << endl;

    for (auto x: P) {
      auto I = kdTreeRange(P, tree, x, k);
      // cout << I.size() << endl;
      // for (auto x: I) cout << x << endl;

      // auto I2 = bruteforceRange<dim>(P, P[0], 4);
      // cout << endl;
      // for (auto x: I2) cout << x << endl;
    }

    cout << "query-time = " << t.get_next() << endl;
  }
  t.stop();
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  double k = P.getOptionIntValue("-k", 1.0);
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
