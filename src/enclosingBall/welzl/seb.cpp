#include <vector>
#include "parlay/sequence.h"
// #include "pbbs/gettime.h"
// #include "pbbs/utils.h"
// #include "pbbs/randPerm.h"
// #include "prefix.h"
// #include "miniDisc.h"
// #include "geometry.h"
#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/welzl/welzl.h"
#include "enclosingBall/ball.h"
// #include "check.h"

//todo make point generic
template<int dim>
pargeo::seb::ball<dim>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P) {
// }
// template<int dim>
// void miniDisc(point<dim>* P, intT n) {
  typedef pargeo::point<dim> pointT;
  typedef pargeo::seb::ball<dim> ballT;

  // static const bool preprocess = false;
  // static const bool serial = false;
  /*
    - 0: plain // do this one for now (todo)
    - 1: mtf
    - 2: pivot + mtf
   */
  // static const int method = 2;

// #ifndef SILENT
//   cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;
// #endif

  // timing t0;t0.start();
//   if(preprocess) {
//     randPerm(P, n);
// #ifndef SILENT
//     cout << "preprocess-time = " << t0.next() << endl;
// #endif
//   }
  //auto support = std::vector<pointT>();
  auto support = parlay::sequence<pointT>();
  ballT D = miniDiscPlain(P, support, ballT());
  return D;

//   ballT D;
//   switch (method) {
//   case 0: {
// #ifndef SILENT
//     cout << "method = plain" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscPlainSerial(P, n, support, ballT());
//     else
//       D = miniDiscPlain(P, n, support, ballT());
//     break;
//   }
//   case 1: {
// #ifndef SILENT
//     cout << "method = mtf" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscMtfSerial(P, n, support, ballT());
//     else
//       D = miniDiscMtf(P, n, support, ballT());
//     break;
//   }
//   case 2: {
// #ifndef SILENT
//     cout << "method = mtf+pivot" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscPivotSerial(P, n, support, ballT());
//     else
//       D = miniDiscPivot(P, n, support, ballT());
//     break;
//   }
//   default:
//     cout << "invalid method, abort" << endl; abort();
//   }

// #ifndef SILENT
//   cout << "seb-time = " << t0.stop() << endl;
//   cout << D.radius() << ", center = " << D.center() << endl;
//   cout << endl;
//   check<dim,ballT>(&D, P, n);
// #else
//   cout << t0.stop() << endl;
// #endif
}


// using namespace std;

// template<int dim>
// ball<dim> miniDiscMtfSerial(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
//   typedef ball<dim> ballT;
//   typedef point<dim> pointT;

//   B = support2Ball(P, support);

//   if (B.size() == dim+1) {
//     return B;
//   }

//   for (intT i=0; i<n; ++i) {

//     if (!B.contain(P[i])) {
//       if (support.size() == B.size()) B.grow(P[i]);
//       else B = ballT();
//       support.push_back(P[i]);
//       B = miniDiscMtfSerial(P, i, support, B);
//       support.pop_back();

//       if (i > dim-support.size()) {
//         swap(P[dim-support.size()], P[i]);}
//     }
//   }

//   return B;
// }

// template<int dim>
// ball<dim> miniDiscMtf(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B, intT* flag=NULL) {
//   typedef ball<dim> ballT;
//   typedef point<dim> pointT;

//   if (n < 2000) return miniDiscMtfSerial(P, n, support, B);

//   B = support2Ball(P, support);
//   if (B.size() == dim+1) return B;

//   bool freeFlag = false;
//   if (!flag) {
//     freeFlag = true;
//     flag = newA(intT, n+1);}

//   auto process = [&](pointT p) {
//                    if (!B.contain(p)) return true;
//                    else return false;
//                  };
//   auto cleanUp = [&](pointT* A, intT i) {
//                    if (support.size() == B.size()) B.grow(A[i]);
//                    else B = ballT();
//                    support.push_back(A[i]);
//                    B = miniDiscMtf(A, i, support, B, flag);
//                    support.pop_back();

//                    if (i > dim-support.size()) {
//                      swap(P[dim-support.size()], P[i]);}
//                  };
//   parallel_prefix(P, n, process, cleanUp, freeFlag, flag);

//   if(freeFlag) free(flag);
//   return B;
// }

// template<int dim>
// intT findPivot(point<dim>* P, intT n, ball<dim> B, intT s) {
//   floatT rSqr = B.radius() * B.radius();
//   floatT dMax = 0;
//   intT bestI = -1;
//   for (intT ii=s; ii<n; ++ii) {
//     floatT tmp = P[ii].distSqr(B.center());
//     if (tmp - rSqr > dMax) { // ||p-c||^2 - r^2
//       bestI = ii;
//       dMax = tmp - rSqr;
//     }
//   }
//   return bestI;
// }

// template<int dim>
// ball<dim> miniDiscPivotSerial(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
//   typedef ball<dim> ballT;
//   typedef point<dim> pointT;

//   B = support2Ball(P, support);

//   if (B.size() == dim+1) {
//     return B;
//   }

//   for (intT i=0; i<n; ++i) {
//     if (!B.contain(P[i])) {
//       intT ii = findPivot(P, n, B, i+1);
//       if (ii > i) swap(P[ii], P[i]);

//       if (support.size() == B.size()) B.grow(P[i]);
//       else B = ballT();
//       support.push_back(P[i]);
//       B = miniDiscMtfSerial(P, i, support, B);
//       support.pop_back();

//       if (i > dim-support.size()) {
//         swap(P[dim-support.size()], P[i]);}
//     }
//   }

//   return B;
// }

// template<int dim>
// ball<dim> miniDiscPivot(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B, intT* flag=NULL) {
//   typedef ball<dim> ballT;
//   typedef point<dim> pointT;

//   if (n < 2000) return miniDiscPivotSerial(P, n, support, B);

//   B = support2Ball(P, support);
//   if (B.size() == dim+1) return B;

//   bool freeFlag = false;
//   if (!flag) {
//     freeFlag = true;
//     flag = newA(intT, n+1);}

//   auto process = [&](pointT p) {
//                    if (!B.contain(p)) return true;
//                    else return false;
//                  };
//   auto cleanUp = [&](pointT* A, intT i) {
//                    intT ii = findPivot(A, i, B, i+1);
//                    if (ii > i) swap(P[ii], P[i]);

//                    if (support.size() == B.size()) B.grow(A[i]);
//                    else B = ballT();
//                    support.push_back(A[i]);
//                    B = miniDiscMtf(A, i, support, B, flag);
//                    support.pop_back();

//                    if (i > dim-support.size()) {
//                      swap(P[dim-support.size()], P[i]);}
//                  };
//   parallel_prefix(P, n, process, cleanUp, freeFlag, flag);

//   if(freeFlag) free(flag);
//   return B;
// }

// template<int dim>
// void miniDisc(point<dim>* P, intT n) {
//   typedef point<dim> pointT;
//   typedef ball<dim> ballT;

//   static const bool preprocess = false;
//   static const bool serial = false;
//   /*
//     - 0: plain
//     - 1: mtf
//     - 2: pivot + mtf
//    */
//   static const int method = 2;

// #ifndef SILENT
//   cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;
// #endif

//   timing t0;t0.start();
//   if(preprocess) {
//     randPerm(P, n);
// #ifndef SILENT
//     cout << "preprocess-time = " << t0.next() << endl;
// #endif
//   }

//   ballT D;
//   switch (method) {
//   case 0: {
// #ifndef SILENT
//     cout << "method = plain" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscPlainSerial(P, n, support, ballT());
//     else
//       D = miniDiscPlain(P, n, support, ballT());
//     break;
//   }
//   case 1: {
// #ifndef SILENT
//     cout << "method = mtf" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscMtfSerial(P, n, support, ballT());
//     else
//       D = miniDiscMtf(P, n, support, ballT());
//     break;
//   }
//   case 2: {
// #ifndef SILENT
//     cout << "method = mtf+pivot" << endl;
// #endif
//     auto support = vector<pointT>();
//     if (serial)
//       D = miniDiscPivotSerial(P, n, support, ballT());
//     else
//       D = miniDiscPivot(P, n, support, ballT());
//     break;
//   }
//   default:
//     cout << "invalid method, abort" << endl; abort();
//   }

// #ifndef SILENT
//   cout << "seb-time = " << t0.stop() << endl;
//   cout << D.radius() << ", center = " << D.center() << endl;
//   cout << endl;
//   check<dim,ballT>(&D, P, n);
// #else
//   cout << t0.stop() << endl;
// #endif
// }

template
pargeo::seb::ball<2>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*> P);

template
pargeo::seb::ball<3>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*> P);

template
pargeo::seb::ball<4>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<4>*, pargeo::point<4>*> P);

template
pargeo::seb::ball<5>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<5>*, pargeo::point<5>*> P);

template
pargeo::seb::ball<6>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<6>*, pargeo::point<6>*> P);

template
pargeo::seb::ball<7>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<7>*, pargeo::point<7>*> P);
