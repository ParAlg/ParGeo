#include <iostream>
#include <algorithm>

#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/getTime.h"
#include "pargeo/parseCommandLine.h"

#include "closestPair/closestPair.h"

template<int dim>
void callClosestPair(parlay::sequence<pargeo::point<dim>>& P) {

  // Perturb points to avoid duplicates
  parlay::parallel_for(0, P.size(), [&](size_t i) {
      for (int j = 0; j < dim; ++ j) {
	double myRand = P[i][j] / 10000;
	P[i][j] += -myRand + 2 * myRand * parlay::hash64(i) /
	  std::numeric_limits<size_t>::max();
      }});

  pargeo::timer t0; t0.start();

  auto cpair = pargeo::closestPair<dim>(P);

  std::cout << "closest pair = " << cpair.first << ", " << cpair.second << "\n";

  std::cout << "time = " << t0.stop() << "\n";
}

int main(int argc, char* argv[]) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  commandLine P(argc, argv, "<inFile>");
  char* iFile = P.getArgument(0);

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> P =
      readPointsFromFile<pargeo::point<2>>(iFile);
    callClosestPair<2>(P);}
  else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> P =
      readPointsFromFile<pargeo::point<3>>(iFile);
    callClosestPair<3>(P);}
  else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> P =
      readPointsFromFile<pargeo::point<4>>(iFile);
    callClosestPair<4>(P);}
  else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> P =
      readPointsFromFile<pargeo::point<5>>(iFile);
    callClosestPair<5>(P);}
  else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> P =
      readPointsFromFile<pargeo::point<6>>(iFile);
    callClosestPair<6>(P);}
  else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> P =
      readPointsFromFile<pargeo::point<7>>(iFile);
    callClosestPair<7>(P);}
  else {
    throw std::runtime_error("Dimension not supported.");
  }
}
