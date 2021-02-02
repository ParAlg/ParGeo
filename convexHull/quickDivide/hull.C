#include "quick.h"
#include "quickDivide.h"

_seq<intT> hull(point2d* P, intT n) {
  timing t; t.start();
  //auto CH = quickHullParallel(P,n);
  auto CH = quickDivide(P,n);
#ifdef SILENT
  cout << t.stop() << endl;
#else
  cout << "hull-time = " << t.stop() << endl;
  check(P, n, CH.A, CH.n);
#endif
  return CH;
}
