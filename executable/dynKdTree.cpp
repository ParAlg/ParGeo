#include <algorithm>
#include <assert.h>
#include <iostream>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "pargeo/dynKdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;
using namespace pargeo::IO;

template<int dim>
void timeKnn(parlay::sequence<pargeo::point<dim>> &P, size_t k, char const *outFile) {
  using namespace pargeo::dynKdTree;

  int n = P.size();

  std::cout << "Test dynamic kdTree of dim " << dim << " on ";
  std::cout << "data of size " << n << "\n";

  pargeo::timer t;
  t.start();

  unique_ptr<rootNode<dim, point<dim>>>
    tree(new rootNode<dim, point<dim>>(P, 0, n / 2));

  std::cout << "build-time = " << t.get_next() << "\n";

  tree->insert(P, n / 2, n);

  std::cout << "size-after-insert = " << tree->size() << "\n";

  std::cout << "insert-time = " << t.get_next() << "\n";

  parlay::parallel_for (0, n, [&](size_t i) {
    tree->kNN(P[i], 5);
  });

  std::cout << "knn-time = " << t.get_next() << "\n";

  tree->erase(P, n / 2, n);

  std::cout << "size-after-erase = " << tree->size() << "\n";

  std::cout << "erase-time = " << t.get_next() << "\n";

  parlay::parallel_for (0, n, [&](size_t i) {
    tree->kNN(P[i], 5);
  });

  std::cout << "knn-time = " << t.get_next() << "\n";

  assert(tree->check());
  assert(tree->size() == n / 2);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-k <param>] [-o <outFile>] <inFile>");
  char* iFile = P.getArgument(0);
  size_t k = P.getOptionIntValue("-k",1);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> Points = readPointsFromFile<pargeo::point<2>>(iFile);
    timeKnn<2>(Points, k, oFile);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> Points = readPointsFromFile<pargeo::point<3>>(iFile);
    timeKnn<3>(Points, k, oFile);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> Points = readPointsFromFile<pargeo::point<4>>(iFile);
    timeKnn<4>(Points, k, oFile);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> Points = readPointsFromFile<pargeo::point<5>>(iFile);
    timeKnn<5>(Points, k, oFile);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> Points = readPointsFromFile<pargeo::point<6>>(iFile);
    timeKnn<6>(Points, k, oFile);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> Points = readPointsFromFile<pargeo::point<7>>(iFile);
    timeKnn<7>(Points, k, oFile);
  } else {
    cout << "dim = " << dim << endl;
    throw std::runtime_error("dimension not supported");
  }
}
