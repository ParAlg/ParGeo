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
#include "geometry.h"
#include "geometryIO.h"
#include "delaunay.h"
#include "pbbs/parallel.h"
#include "pbbs/gettime.h"
#include "pbbs/parseCommandLine.h"
using namespace std;
using namespace benchIO;

// *************************************************************
//  TIMING
// *************************************************************

uint64_t hash64(uint64_t u)
{
  uint64_t v = u * 3935559000370003845 + 2691343689449507681;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >>  4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v <<  5;
  return v;
}

template<int dim>
void timeDelaunay(point2d* pts, intT n, int rounds, char* outFile, int perturb) {
  cout << "nmax = " << n << endl;
  if (perturb) {
    par_for(intT i=0; i<n; ++i) {
      double myRand = pts[i][0] / 10000;
      pts[i].p.x[0] += -myRand + 2*myRand*hash64(i)/numeric_limits<unsigned long>::max();
      myRand = pts[i][1] / 10000;
      pts[i].p.x[1] += -myRand + 2*myRand*hash64(i)/numeric_limits<unsigned long>::max();
    }
  }

  triangles<point2d> R;
  for (int i=0; i < rounds; i++) {
    timing t0; t0.start();
    R = delaunay(pts, n);
    cout << "timing = " << t0.stop() << endl;
  }
  cout << endl;
  if (outFile != NULL) {
    cout << "writing triangles to file" << endl;
    writeTrianglesToFile(R, outFile);
  }
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-p <perturb points? 0/1>] [-csv #columns> [-nmax <#points>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int perturb = P.getOptionIntValue("-p",0);
  int nMax = P.getOptionIntValue("-nmax",-1);
  int csvCol = P.getOptionIntValue("-csv",-1);
  bool readCsv = csvCol > 0;

  int dim;
  if(!readCsv) dim = readPointsDimensionFromFile(iFile);
  else dim = csvCol;

  printScheduler();
  cout << "perturb points = " << perturb << endl;

  if (dim == 1) {
    cout << "dimension 1 not supported, abort" << endl; abort();
  } else if (dim == 2) {
    _seq<point<2>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<2>>(iFile);
    else PIn = readPointsFromFileCSV<point<2>>(iFile, csvCol);
    point2d* PP = newA(point2d, PIn.n);
    par_for(intT i=0; i<PIn.n; ++i) PP[i] = point2d(PIn.A[i]);//convert to 2d points
    timeDelaunay<2>(PP, nMax>0? nMax : PIn.n, rounds, oFile, perturb);
    free(PP);
  }
  else {
    cout << "dimension " << dim << " not yet supported, abort" << endl; abort();
  }
}
