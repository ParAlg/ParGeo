#include <iostream>
#include <algorithm>
#include "spanner.h"
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
void timeSpanner(point<dim>* pts, intT n, int t, int rounds, char* outFile, int perturb) {
  cout << "nmax = " << n << endl;
  if (perturb) {
    parallel_for(0, n, [&](intT i) {
	for (int j = 0; j < dim; ++ j) {
	  double myRand = pts[i][j] / 10000;
	  pts[i].x[j] += -myRand + 2*myRand*hash64(i)/numeric_limits<unsigned long>::max();
	}
      });
  }

  for (int i=0; i < rounds; i++) {
    timing t0; t0.start();
    auto R = spanner<dim>(pts, n, t);
    cout << "timing = " << t0.stop() << endl;
    floatT sum = 0;
    for (int j=0; j<n-1; j++) {
      sum += R[j].weight;
    }
    cout << "weight-sum = " << sum << endl;
  }
  cout << endl;

  // if (outFile != NULL) write(R); // todo
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] [-r <rounds>] [-p <perturb points? 0/1>] [-csv <#columns> [-nmax <#points>] <inFile>");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int rounds = P.getOptionIntValue("-r",1);
  int perturb = P.getOptionIntValue("-p",0);
  int t = P.getOptionIntValue("-t",2);//t-spanner parameter
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
    timeSpanner<2>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 3) {
    _seq<point<3>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<3>>(iFile);
    else PIn = readPointsFromFileCSV<point<3>>(iFile, csvCol);
    timeSpanner<3>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 4) {
    _seq<point<4>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<4>>(iFile);
    else PIn = readPointsFromFileCSV<point<4>>(iFile, csvCol);
    timeSpanner<4>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 5) {
    _seq<point<5>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<5>>(iFile);
    else PIn = readPointsFromFileCSV<point<5>>(iFile, csvCol);
    timeSpanner<5>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 6) {
    _seq<point<6>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<6>>(iFile);
    else PIn = readPointsFromFileCSV<point<6>>(iFile, csvCol);
    timeSpanner<6>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 7) {
    _seq<point<7>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<7>>(iFile);
    else PIn = readPointsFromFileCSV<point<7>>(iFile, csvCol);
    timeSpanner<7>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 8) {
    _seq<point<8>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<8>>(iFile);
    else PIn = readPointsFromFileCSV<point<8>>(iFile, csvCol);
    timeSpanner<8>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else if (dim == 9) {
    _seq<point<9>> PIn;
    if(!readCsv) PIn = readPointsFromFile<point<9>>(iFile);
    else PIn = readPointsFromFileCSV<point<9>>(iFile, csvCol);
    timeSpanner<9>(PIn.A, nMax>0? nMax : PIn.n, t, rounds, oFile, perturb);
  }
  else {
    cout << "dimension " << dim << " not yet supported, abort" << endl; abort();
  }
}
