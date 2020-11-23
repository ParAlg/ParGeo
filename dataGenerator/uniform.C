#include <math.h>
#include <string>
#include "pbbs/parallel.h"
#include "pbbs/parseCommandLine.h"
#include "pbbs/randPerm.h"
#include "IO.h"
#include "geometry.h"
#include "geometryIO.h"
#include "geometryData.h"
#include "dataGen.h"
using namespace benchIO;
using namespace dataGen;
using namespace std;

int main(int argc, char* argv[]) {
  string text = "[-s] [-S] [-p] [-d {2--9}] [-csv] n <outFile>";
  text += "\n -s: in sphere";
  text += "\n -S: on sphere";
  text += "\n -p: on parabola/parabloid (2d/3d only)";
  text += "\n default : in cube";
  commandLine P(argc,argv,text);
  pair<intT,char*> in = P.sizeAndFileName();
  intT n = in.first;
  char* fname = in.second;
  int dims = P.getOptionIntValue("-d", 2);
  bool inSphere = P.getOption("-s");
  bool onSphere = P.getOption("-S");
  bool parabloid = P.getOption("-p");
  bool csvOut = P.getOption("-csv");

  if (parabloid && dims > 3) {
    cout << "error, cannot generate >3d parabloid, abort" << endl;
    abort();
  }

  if (dims == 2) {
    point<2>* P = newA(point<2>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<2>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<2>(i, n);
      else if (parabloid) P[i] = randOnParabola(i, n);
      else P[i] = randNd<2>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 3) {
    point<3>* P = newA(point<3>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<3>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<3>(i, n);
      else if (parabloid) P[i] = randOnParabloid(i, n);
      else P[i] = randNd<3>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  else if (dims == 4) {
    point<4>* P = newA(point<4>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<4>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<4>(i, n);
      else P[i] = randNd<4>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 5) {
    point<5>* P = newA(point<5>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<5>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<5>(i, n);
      else P[i] = randNd<5>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 6) {
    point<6>* P = newA(point<6>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<6>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<6>(i, n);
      else P[i] = randNd<6>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 7) {
    point<7>* P = newA(point<7>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<7>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<7>(i, n);
      else P[i] = randNd<7>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 8) {
    point<8>* P = newA(point<8>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<8>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<8>(i, n);
      else P[i] = randNd<8>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  if (dims == 9) {
    point<9>* P = newA(point<9>, n);
    par_for (intT i=0; i<n; i++)
      if (inSphere) P[i] = randInUnitSphere<9>(i, sqrt(float(n)));
      else if (onSphere) P[i] = randOnUnitSphere<9>(i, n);
      else P[i] = randNd<9>(i, sqrt(float(n)));
    randPerm(P, n);
    if(csvOut) return writePointsToFileCSV(P, n, fname);
    else return writePointsToFile(P, n, fname);
  }
  cout << "bad argument" << endl;
  P.badArgument();
  return 1;
}
