#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "common/get_time.h"
#include "common/geometry.h"
#include "common/geometryIO.h"
#include "common/parse_command_line.h"
#include "geometry/point.h"
#include "hull.h"

using namespace std;
using namespace benchIO;

template <class pt>
void timeHull(parlay::sequence<pt> &P, size_t numProc, int rounds, char const *outFile) {
  timer t; t.start();
  for(int i=0; i<rounds; ++i) {
    hull3d(P, numProc);
    cout << "round-time = " << t.get_next() << endl;
  }
  t.stop();
  //if (outFile != NULL) writeIntSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int numProc = P.getOptionIntValue("-p",1);
  int rounds = P.getOptionIntValue("-r",1);

  int dim = readDimensionFromFile(iFile);//todo make cheaper
  if (dim != 3) {
    cout << "Error, convexHull3D only takes 3d inputs, abort." << endl;
    abort();
  }

  parlay::sequence<pargeo::fpoint<3>> Points = readPointsFromFile<pargeo::fpoint<3>>(iFile);
  timeHull<pargeo::fpoint<3>>(Points, numProc, rounds, oFile);
}
