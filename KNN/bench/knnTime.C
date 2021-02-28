#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "common/get_time.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "knn.h"

using namespace std;
using namespace benchIO;

using coord = double;

template<int dim>
void timeKnn(parlay::sequence<point<dim>> &P, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    knn<dim>(P);
    cout << "round-time = " << t.get_next() << endl;
  }
  t.stop();
  //if (outFile != NULL) writeIntSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  int dim = readDimensionFromFile(iFile);//todo make cheaper

  if (dim == 2) {
    parlay::sequence<point<2>> Points = readPointsFromFile<point<2>>(iFile);
    timeKnn<2>(Points, rounds, oFile);
  } else if (dim == 3) {
    parlay::sequence<point<3>> Points = readPointsFromFile<point<3>>(iFile);
    timeKnn<3>(Points, rounds, oFile);
  }
}
