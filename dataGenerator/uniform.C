#include <math.h>
#include "pbbs/parallel.h"
#include "IO.h"
#include "geometry.h"
#include "geometryIO.h"
#include "geometryData.h"
#include "dataGen.h"
#include "pbbs/parseCommandLine.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-csv] [-s] [-S] [-d {2--9}] n <outFile>");
  pair<intT,char*> in = P.sizeAndFileName();
  intT n = in.first;
  char* fname = in.second;
  int dims = P.getOptionIntValue("-d", 2);
  bool inSphere = P.getOption("-s");
  bool onSphere = P.getOption("-S");
  bool csvOut = P.getOption("-csv");

  if (dims == 2) {
    point<2>* P = newA(point<2>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<2>(i);
      else if (onSphere) P[i] = randOnUnitSphere<2>(i);
      else P[i] = randNd<2>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 3) {
    point<3>* P = newA(point<3>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<3>(i);
      else if (onSphere) P[i] = randOnUnitSphere<3>(i);
      else P[i] = randNd<3>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  else if (dims == 4) {
    point<4>* P = newA(point<4>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<4>(i);
      else if (onSphere) P[i] = randOnUnitSphere<4>(i);
      else P[i] = randNd<4>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 5) {
    point<5>* P = newA(point<5>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<5>(i);
      else if (onSphere) P[i] = randOnUnitSphere<5>(i);
      else P[i] = randNd<5>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 6) {
    point<6>* P = newA(point<6>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<6>(i);
      else if (onSphere) P[i] = randOnUnitSphere<6>(i);
      else P[i] = randNd<6>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 7) {
    point<7>* P = newA(point<7>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<7>(i);
      else if (onSphere) P[i] = randOnUnitSphere<7>(i);
      else P[i] = randNd<7>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 8) {
    point<8>* P = newA(point<8>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<8>(i);
      else if (onSphere) P[i] = randOnUnitSphere<8>(i);
      else P[i] = randNd<8>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 9) {
    point<9>* P = newA(point<9>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<9>(i);
      else if (onSphere) P[i] = randOnUnitSphere<9>(i);
      else P[i] = randNd<9>(i);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  cout << "bad argument" << endl;
  P.badArgument();
  return 1;
}
