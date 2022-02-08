#include <iostream>
#include <algorithm>
#include "euclideanMst/euclideanMst.h"
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/graphIO.h"
#include "pargeo/parseCommandLine.h"
#include "spatialGraph/spatialGraph.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::pointIO;
using namespace pargeo::graphIO;

template<int dim>
void callAlgo(int algo, parlay::sequence<pargeo::point<dim>> &P, double param, char const *outFile) {

  auto dimLimit2 = [&]() {
    if (dim != 2)
      throw std::runtime_error("Error, only dimension 2 is supported for this generator.");};

  auto dimLimit7 = [&]() {
    if (dim < 2 || dim > 7)
      throw std::runtime_error("Error, only dimension 2-7 is supported for this generator.");};

  switch(algo) {
  case 1: // kNN graph
    {
      dimLimit7();
      auto I = knnGraph<dim>(P, size_t(param));
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  case 2: // Delaunay graph
    {
      dimLimit2();
      auto I = delaunayGraph<dim>(P);
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  case 3: // Beta skeleton
    {
      dimLimit2();
      auto I = pargeo::betaSkeleton<dim>(P, param);
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  case 4: // Gabriel graph
    {
      dimLimit2();
      auto I = gabrielGraph<dim>(P);
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  case 5: // Euclidean MST
    {
      dimLimit7();
      parallel_for(0, P.size(), [&](size_t i) { // point perturbation
	for (int j = 0; j < dim; ++ j) {
	  double myRand = P[i][j] / 10000;
	  P[i][j] += -myRand + 2*myRand*hash64(i)/numeric_limits<size_t>::max();
	}});
      auto I = euclideanMst<dim>(P);
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  case 6: // WSPD spanner
    {
      dimLimit7();
      parallel_for(0, P.size(), [&](size_t i) { // point perturbation
	for (int j = 0; j < dim; ++ j) {
	  double myRand = P[i][j] / 10000;
	  P[i][j] += -myRand + 2*myRand*hash64(i)/numeric_limits<size_t>::max();
	}});
      auto I = spanner<dim>(P, param);
      if (outFile != NULL) graphIO::writeEdgeSeqToFile(I, outFile);
      break;
    }
  default:
    {
      throw std::runtime_error("Error, invalid input generator.");
    }
  }

}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-algo <algorithm>] [-param <param>] [-o <outFile>] <inFile>\n\nSpatial graph generators:\n\n kNN Graph (2-7d):     -algo 1 -param <k>\n\n Delaunay Graph (2d):  -algo 2\n\n Beta Skeleton (2d):   -algo 3 -param <beta>\n\n Gabriel Graph (2d):   -algo 4\n\n Euclidean MST (2-7d): -algo 5\n\n WSPD Spanner (2-7d):  -algo 6 -param <t>\n\nE.g. 1-NN graph   ./graphGenerator -algo 1 -param 1 -o graph.txt input.csv");
  char* iFile = P.getArgument(0);
  int algo = P.getOptionIntValue("-algo",1);
  double param = P.getOptionDoubleValue("-param",1.0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> P = readPointsFromFile<pargeo::point<2>>(iFile);
    callAlgo<2>(algo, P, param, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> P = readPointsFromFile<pargeo::point<3>>(iFile);
    callAlgo<3>(algo, P, param, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> P = readPointsFromFile<pargeo::point<4>>(iFile);
    callAlgo<4>(algo, P, param, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> P = readPointsFromFile<pargeo::point<5>>(iFile);
    callAlgo<5>(algo, P, param, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> P = readPointsFromFile<pargeo::point<6>>(iFile);
    callAlgo<6>(algo, P, param, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> P = readPointsFromFile<pargeo::point<7>>(iFile);
    callAlgo<7>(algo, P, param, oFile);
  } else {
    throw std::runtime_error("Unsupported dimensionality > 7.");
  }
}
