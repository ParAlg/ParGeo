#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "common/get_time.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "hull.h"

using namespace std;
using namespace benchIO;

using coord = double;

void timeHull(parlay::sequence<point<3>> const &P, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    hull(P);
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

  if (dim == 3) {
    parlay::sequence<point<3>> Points = readPointsFromFile<point<3>>(iFile);
    timeHull(Points, rounds, oFile);
  } else {
    cout << "Error, convexHull3D only takes 3d inputs, abort." << endl;
    abort();
  }
}
