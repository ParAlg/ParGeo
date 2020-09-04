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
#include "closestPair.h"
#include "geometry.h"
#include "geometryIO.h"
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
void timeClosestPair(point<dim>* pts, intT n, int batches, int rounds, char* outFile, int perturb) {
  cout << "nmax = " << n << endl;
  if (perturb) {
    par_for(intT i=0; i<n; ++i) {
      for (int j = 0; j < dim; ++ j) {
        double myRand = pts[i][j] / 10000;
        //double myRand = pts[i][j] / 10;
        pts[i].x[j] += -myRand + 2*myRand*hash64(i)/numeric_limits<unsigned long>::max();
      }
    }
  }

  pair<point<dim>, point<dim>> R;
  for (int i=0; i < rounds; i++) {
    timing t0; t0.start();
    R = closestPair<dim>(pts, n, batches);
    cout << "timing = " << t0.stop() << endl;
  }
  cout << endl;

  // if (outFile != NULL) write(R); // todo
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-b <#batches>] [-p <perturb points? 0/1>] [-csv <csv total #columns> -scol <csv start column inclusive> -ecol <csv end column exclusive>] [-nmax <#points to feed>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int perturb = P.getOptionIntValue("-p",0);
  int nMax = P.getOptionIntValue("-nmax",-1);
  int csvCol = P.getOptionIntValue("-csv",-1);
  int csvSCol = P.getOptionIntValue("-scol",-1);
  int csvECol = P.getOptionIntValue("-ecol",-1);
  bool readCsv = csvCol > 0;
  int batches = P.getOptionIntValue("-b",1);

  int dim;
  if(!readCsv) dim = readPointsDimensionFromFile(iFile);
  else dim = csvECol - csvSCol;

  printScheduler();
  cout << "perturb points = " << perturb << endl;

  if (dim == 1) {
    cout << "dimension 1 not supported, abort" << endl; abort();
  } else if (dim == 2) {
    _seq<point<2>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<2>>(iFile);
    else PIn = readPointsFromFileCSV<point<2>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<2>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 3) {
    _seq<point<3>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<3>>(iFile);
    else PIn = readPointsFromFileCSV<point<3>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<3>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 4) {
    _seq<point<4>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<4>>(iFile);
    else PIn = readPointsFromFileCSV<point<4>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<4>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 5) {
    _seq<point<5>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<5>>(iFile);
    else PIn = readPointsFromFileCSV<point<5>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<5>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 6) {
    _seq<point<6>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<6>>(iFile);
    else PIn = readPointsFromFileCSV<point<6>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<6>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 7) {
    _seq<point<7>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<7>>(iFile);
    else PIn = readPointsFromFileCSV<point<7>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<7>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 8) {
    _seq<point<8>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<8>>(iFile);
    else PIn = readPointsFromFileCSV<point<8>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<8>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else if (dim == 9) {
    _seq<point<9>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<9>>(iFile);
    else PIn = readPointsFromFileCSV<point<9>>(iFile, csvCol, csvSCol, csvECol);
    timeClosestPair<9>(PIn.A, nMax>0? nMax : PIn.n, batches, rounds, oFile, perturb);
  }
  else {
    cout << "dimension " << dim << " not yet supported, abort" << endl; abort();
  }
}
