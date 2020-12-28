#ifndef HULL_H
#define HULL_H

#include "geometry.h"
#include "pbbs/sequence.h"

_seq<intT> hull(point2d*, intT);

static const floatT numericKnob = 1e-9;

struct aboveLine2 {
  point2d l, r;
  point2d* P;
  aboveLine2(point2d* _P, point2d _l, point2d _r) : P(_P), l(_l), r(_r) {}
  bool operator() (intT i) {return triArea(l, r, P[i]) > numericKnob;}
};

inline void check(point2d* P, intT n, intT* I, intT m, point2d* CH=NULL) {
  if (!CH) CH = P;
  // for(intT i=0; i<m; ++i) {
  //   cout << CH[I[i]] << " ";
  // }
  // cout << endl;
  cout << "hull size = " << m << endl;
  intT nume = 0;
  for(intT i=0; i<m; ++i) {
    auto al = aboveLine2(P, CH[I[i]], CH[I[(i+1)%m]]);
    for(intT j=0; j<n; ++j) {
      if(al(j)) {
	floatT error = triArea(CH[I[i]], CH[I[(i+1)%m]], P[j]);
	if (error > 1e-7) {
	  cout << "incorrect hull, point out at " << j << " = " << P[j] << endl;
	  cout << "error = " << error << endl;
	  if(nume>0) cout << "warning, also found " << nume << " small numerical errors" << endl;
	  abort();
	} else nume++;
      }}
  }
  if (nume > 0)
    cout << "warning, found " << nume << " small numerical errors" << endl;
  else
    cout << "hull correct" << endl;
}

#endif
