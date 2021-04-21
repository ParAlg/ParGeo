// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include "pargeo/geometryIO.h"
#include "pargeo/parseCommandLine.h"
#include "delaunayTriangulation/geometry.h"
#include "delaunayTriangulation/triangleIO.h"
#include "delaunayTriangulation/delaunay.h"
#include "parlay/primitives.h"

using namespace std;
using namespace benchIO;
using namespace pbbsbench;

// *************************************************************
//  TIMING
// *************************************************************

void timeDelaunay(parlay::sequence<pointT> &pts, int rounds, char* outFile) {
  triangles<pointT> R;
  for(int i=0; i<rounds; ++i) {
    R = delaunay(pts);
  }
  cout << endl;
  if (outFile != NULL) writeTrianglesToFile(R, outFile);
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);

  parlay::sequence<pointT> PI = readPointsFromFile<pointT>(iFile);
  timeDelaunay(PI, rounds, oFile);
}
