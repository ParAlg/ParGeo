#include <iostream>
#include <algorithm>
#include "parlay/utilities.h"
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/getTime.h"
#include "pargeo/parseCommandLine.h"
#include "clustering/hdbscan.h"
#include "clustering/dendrogram.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::pointIO;

// *************************************************************
//  TIMING
// *************************************************************

template<int dim>
void timeHdbscan(sequence<point<dim>>& P, size_t minPts, int rounds, char* outFile, int perturb) {
  if (perturb) {
    parallel_for(0, P.size(), [&](size_t i) {
	for (int j = 0; j < dim; ++ j) {
	  double myRand = P[i][j] / 10000;
	  P[i][j] += -myRand + 2*myRand*hash64(i)/numeric_limits<size_t>::max();
	}});
  }

  for (int i=0; i < rounds; i++) {
    timer t0; t0.start();
    sequence<wghEdge> E = hdbscan<dim>(P, minPts);
    dendrogram(E, P.size());
    cout << "timing = " << t0.stop() << endl;
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-p <0/1 perturb points>] <inFile>");
  char* iFile = P.getArgument(0);
  size_t k = P.getOptionIntValue("-k",1);
  size_t minPts = P.getOptionIntValue("-m",1);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int perturb = P.getOptionIntValue("-p",0);

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeHdbscan<2>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeHdbscan<3>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeHdbscan<4>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeHdbscan<5>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeHdbscan<6>(Points, minPts, rounds, oFile, perturb);}
  else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeHdbscan<7>(Points, minPts, rounds, oFile, perturb);}
  else {
    cout << "dimension " << dim << " not yet supported, abort" << endl; abort();
  }
}
