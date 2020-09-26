#ifndef CHECK_H
#define CHECK_H

#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "geometry.h"

// *************************************************************
//    CHECKER
// *************************************************************

template<int dim>
bool check(sphere<dim>* circle, point<dim>* P, intT n, bool verbose=true) {
  //code for verifying correctness
  for (intT i=0; i<n; ++i) {
    if(!circle->contain(P[i])) {
      if (verbose) {
        cout << "outside point = " << P[i] << endl;
        cout << "dist = " << P[i].pointDist(circle->center()) << endl;
        abort();
      } else {
        return false;
      }
    }
  }
  if(verbose) cout << "correctness verified" << endl;
  return true;
}

#endif
