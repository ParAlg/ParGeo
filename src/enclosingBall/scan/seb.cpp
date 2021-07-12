#include <vector>
#include "parlay/sequence.h"
#include "enclosingBall/scan/seb.h"
#include "enclosingBall/welzl/welzl.h"
#include "enclosingBall/ball.h"

//todo make point generic

template<int dim>
size_t quadrant(pargeo::point<dim> p, pargeo::point<dim> center) {
  int index = 0;
  int offset = 1;
  for (int i = 0; i < dim; ++ i) {
    if (p[i] > center[i]) index += offset;
    offset *= 2;
  }
  return index;
}

template<int dim>
bool ortScanSerial(pargeo::point<dim> c,
		   typename pargeo::point<dim>::floatT rSqr,
		   parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
		   parlay::sequence<pargeo::point<dim>>& support,
		   typename pargeo::point<dim>::floatT* dist) {

  using pointT = pargeo::point<dim>;
  using floatT = typename pargeo::point<dim>::floatT;

  int dd = int(pow(2.0, dim));
  long idx[dd];
  for(int i=0; i<dd; ++i) idx[i] = -1;

  for (size_t i = 0; i < P.size(); ++ i) {
    floatT dSqr = P[i].distSqr(c);
    if (dSqr > rSqr + c.eps) {
      int o = quadrant(c, P[i]);
      if (dSqr > dist[o] + c.eps) {
        dist[o] = dSqr;
        idx[o] = i;
      }
    }
  }

  bool hasOut = false;
  for(int i = 0; i < dd; ++ i) {
    if (idx[i] != -1) {
      hasOut = true;
      support.push_back(P[idx[i]]);
    }
  }
  return hasOut;
}

template<int dim>
bool ortScan(pargeo::point<dim> c,
	     typename pargeo::point<dim>::floatT rSqr,
	     parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> A,
	     parlay::sequence<pargeo::point<dim>>& support,
	     typename pargeo::point<dim>::floatT* distGlobal) {
  using pointT = pargeo::point<dim>;
  using floatT = typename pargeo::point<dim>::floatT;

  if (A.size() < 2000) return ortScanSerial(c, rSqr, A, support, distGlobal);

  int dd = int(pow(2.0, dim));
  int P = std::max((size_t)96, parlay::num_workers());
  size_t blockSize = (A.size() + P - 1) / P;
  long idx[dd * P];
  floatT dist[dd * P];
  for (size_t i = 0; i < dd * P; ++ i) idx[i] = -1;
  for (int i = 0; i < P; ++ i) {
    for (int j = 0; j < dd; ++ j) {
      dist[i * dd + j] = distGlobal[j];}
  }

  parlay::parallel_for(0, P,
	       [&](int p) {
		 size_t s = p * blockSize;
		 size_t e = std::min((size_t)(p + 1) * blockSize, A.size());
		 long* locIdx = idx + p*dd;
		 floatT* locDist = dist + p * dd;
		 for (size_t i = s; i < e; ++ i) {
		   floatT dSqr = A[i].distSqr(c);
		   if (dSqr > rSqr + c.eps) {
		     int o = quadrant(c, A[i]);
		     if (dSqr > locDist[o] + c.eps) {
		       locDist[o] = dSqr;
		       locIdx[o] = i;}
		   }
		 }
	       }, 1);

  long idxGlobal[dd];
  for (int o = 0; o < dd; ++ o) idxGlobal[o] = -1;

  for(int p = 0; p < P; ++ p) {
    for(int o = 0; o < dd; ++ o) {
      long* locIdx = idx + p * dd;
      floatT* locDist = dist + p * dd;
      if (locIdx[o] != -1) {
        if(locDist[o] > distGlobal[o]) {
          idxGlobal[o] = locIdx[o];
          distGlobal[o] = locDist[o];}}
    }
  }

  bool hasOut = false;
  for(int o = 0; o < dd; ++ o) {
    if (idxGlobal[o] != -1) {
      hasOut = true;
      support.push_back(A[idxGlobal[o]]);}
  }
  return hasOut;
}

template<int dim>
pargeo::seb::ball<dim>
miniDiscOrt(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;
  using floatT = typename pointT::floatT;

  size_t sample = dim*3;
  ballT B;
  if (sample > P.size()) {
    parlay::sequence<pointT> support;
    return welzlParallel<dim>(P.cut(0, sample), support, B);
  } else {
    parlay::sequence<pointT> support;
    B = welzlParallel<dim>(P.cut(0, sample), support, B);
  }

  int dd = int(pow(2.0, dim));
  floatT dist[dd];
  for(int i=0; i<dd; ++i) dist[i] = -1;

  while (1) {
    parlay::sequence<pointT> support;
    for(size_t i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScan<dim>(B.center(), B.radius() * B.radius(), P, support, dist);

    if (!found) {
      return B;
    } else {
      auto supportNew = parlay::sequence<pointT>();
      B = welzlParallel<dim>(parlay::make_slice(support), supportNew, ballT());
      for(int i = 0; i < dd; ++ i) dist[i] = -1;
    }
  }
  return B;
}

template<int dim>
pargeo::seb::ball<dim>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P) {
  return miniDiscOrt(P);
}

template
pargeo::seb::ball<2>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*> P);

template
pargeo::seb::ball<3>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*> P);

template
pargeo::seb::ball<4>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<4>*, pargeo::point<4>*> P);

template
pargeo::seb::ball<5>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<5>*, pargeo::point<5>*> P);

template
pargeo::seb::ball<6>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<6>*, pargeo::point<6>*> P);

template
pargeo::seb::ball<7>
pargeo::seb::scan::compute(parlay::slice<pargeo::point<7>*, pargeo::point<7>*> P);
