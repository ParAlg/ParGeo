#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"
#include "convexHull3d/hull.h"

using namespace std;
using namespace pargeo;
using namespace pargeo::pointIO;

template <class pt>
void timeHull(parlay::sequence<pt> &P, int rounds, char const *outFile) {
#ifndef SILENT
  timer t; t.start();
#endif
  for(int i=0; i<rounds; ++i) {
    hull3d(P);
#ifndef SILENT
    cout << "round-time = " << t.get_next() << endl;
#endif
  }
#ifndef SILENT
  t.stop();
#endif
  //if (outFile != NULL) writeIntSeqToFile(I, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  int dim = readHeader(iFile);
  if (dim != 3) {
    cout << "Error, convexHull3D only takes 3d inputs, abort." << endl;
    abort();
  }

  parlay::sequence<pargeo::fpoint<3>> Points = readPointsFromFile<pargeo::fpoint<3>>(iFile);
  timeHull<pargeo::fpoint<3>>(Points, rounds, oFile);
}
