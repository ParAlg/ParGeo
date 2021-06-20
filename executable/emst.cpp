#include <iostream>
#include <algorithm>
#include "parlay/utilities.h"
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/graphIO.h"
#include "pargeo/getTime.h"
#include "pargeo/parseCommandLine.h"
#include "euclideanMst/euclideanMst.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::pointIO;
using namespace pargeo::graphIO;

template<int dim>
void timeEmst(sequence<point<dim>>& P, char* outFile, int perturb) {
  if (perturb) {
    parallel_for(0, P.size(), [&](size_t i) {
	for (int j = 0; j < dim; ++ j) {
	  double myRand = P[i][j] / 10000;
	  P[i][j] += -myRand + 2*myRand*hash64(i)/numeric_limits<size_t>::max();
	}});
  }

  pargeo::timer t0; t0.start();
  auto I = euclideanMst<dim>(P);
  cout << "time = " << t0.stop() << endl;
  if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int perturb = P.getOptionIntValue("-p",0);

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeEmst<2>(Points, oFile, perturb);}
  else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeEmst<3>(Points, oFile, perturb);}
  else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeEmst<4>(Points, oFile, perturb);}
  else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeEmst<5>(Points, oFile, perturb);}
  else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeEmst<6>(Points, oFile, perturb);}
  else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeEmst<7>(Points, oFile, perturb);}
  else {
    throw std::runtime_error("dimension not yet supported");
  }
}
