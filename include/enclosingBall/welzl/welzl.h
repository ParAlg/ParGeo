#pragma once

#include "pargeo/point.h"
#include "parlay/sequence.h"

#include "enclosingBall/ball.h"
#include "enclosingBall/prefixFor.h"

/*------------ Primitives ------------*/

template<int dim>
pargeo::seb::ball<dim> support2Ball(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
				    parlay::sequence<pargeo::point<dim>>& support) {
  using ballT = pargeo::seb::ball<dim>;

  ballT B;
  if (B.isEmpty()) {
    if (support.size() == 0) {
      B = ballT(P.cut(0, 2));
    } else if (support.size() == 1) {
      support.push_back(P[0]);
      B = ballT(parlay::make_slice(support));
      support.pop_back();
    } else { // >=2
      B = ballT(parlay::make_slice(support));
    }
  }
  return B;
}

template<int dim>
size_t findPivot(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
		 pargeo::seb::ball<dim> B,
		 size_t s) {
  using floatT = typename pargeo::point<dim>::floatT;

  floatT rSqr = B.radius() * B.radius();
  floatT dMax = 0;
  size_t bestI = -1;
  for (size_t ii = s; ii < P.size(); ++ ii) {
    floatT tmp = P[ii].distSqr(B.center());
    if (tmp - rSqr > dMax) { // ||p-c||^2 - r^2                                                                           
      bestI = ii;
      dMax = tmp - rSqr;
    }
  }
  return bestI;
}

/*------------ The vanilla Welzl's algorithm ------------*/

template<int dim>
pargeo::seb::ball<dim> welzlSerial(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					   parlay::sequence<pargeo::point<dim>>& support,
					   pargeo::seb::ball<dim> B) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  B = support2Ball(P, support);

  if (B.size() == dim+1) return B;

  for (size_t i=0; i<P.size(); ++i) {

    if (!B.contain(P[i])) { // process a conflict, ow keep going
      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();

      support.push_back(P[i]);
      // compute a ball on the prefix
      B = welzlSerial(P.cut(0, i), support, B);
      support.pop_back();
    }

  }

  return B;
}

template<int dim>
pargeo::seb::ball<dim> welzlParallel(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
				     parlay::sequence<pargeo::point<dim>>& support,
				     pargeo::seb::ball<dim> B,
				     parlay::sequence<size_t>* flag=NULL) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  if (P.size() < 2000) return welzlSerial<dim>(P, support, B);

  B = support2Ball(P, support);
  if (B.size() == dim+1) return B;

  bool freeFlag = false;
  if (!flag) {
    freeFlag = true;
    flag = new parlay::sequence<size_t>(P.size()+1);
  }

  auto process = [&](pointT p) {
                   if (!B.contain(p)) return true;
                   else return false;
                 };

  auto cleanUp = [&](parlay::slice<pointT*, pointT*> A, size_t i) {
                   if (support.size() == B.size()) B.grow(A[i]);
                   else B = ballT();
                   support.push_back(A[i]);
                   B = welzlParallel(A.cut(0, i), support, B, flag);
                   support.pop_back();
                 };

  parallel_prefix(P, process, cleanUp, flag, 200000, 500000);

  if(freeFlag) delete flag;
  return B;

}

/*------------ Welzl + mtf algorithm ------------*/

template<int dim>
pargeo::seb::ball<dim> welzlMtfSerial(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					   parlay::sequence<pargeo::point<dim>>& support,
					   pargeo::seb::ball<dim> B) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  B = support2Ball(P, support);

  if (B.size() == dim+1) return B;

  for (size_t i=0; i<P.size(); ++i) {

    if (!B.contain(P[i])) { // process a conflict, ow keep going

      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();

      support.push_back(P[i]);

      // compute a ball on the prefix
      B = welzlMtfSerial(P.cut(0, i), support, B);
      support.pop_back();

      if (i > dim-support.size()) {
	std::swap(P[dim-support.size()], P[i]);}
    }

  }

  return B;
}

template<int dim>
pargeo::seb::ball<dim> welzlMtfParallel(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					     parlay::sequence<pargeo::point<dim>>& support,
					     pargeo::seb::ball<dim> B,
					     parlay::sequence<size_t>* flag=NULL) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  if (P.size() < 2000) return welzlMtfSerial<dim>(P, support, B);

  B = support2Ball(P, support);
  if (B.size() == dim+1) return B;

  bool freeFlag = false;
  if (!flag) {
    freeFlag = true;
    flag = new parlay::sequence<size_t>(P.size()+1);
  }

  auto process = [&](pointT p) {
                   if (!B.contain(p)) return true;
                   else return false;
                 };

  auto cleanUp = [&](parlay::slice<pointT*, pointT*> A, size_t i) {

                   if (support.size() == B.size()) B.grow(A[i]);
                   else B = ballT();
                   support.push_back(A[i]);
                   B = welzlMtfParallel(A.cut(0, i), support, B, flag);
                   support.pop_back();

		   if (i > dim-support.size()) {
		     std::swap(P[dim-support.size()], P[i]);}
                 };

  parallel_prefix(P, process, cleanUp, flag, 200000, 500000);

  if(freeFlag) delete flag;
  return B;

}

/*------------ Welzl + pivoting + mtf algorithm ------------*/

template<int dim>
pargeo::seb::ball<dim> welzlMtfPivotSerial(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					   parlay::sequence<pargeo::point<dim>>& support,
					   pargeo::seb::ball<dim> B) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  B = support2Ball(P, support);

  if (B.size() == dim+1) return B;

  for (size_t i=0; i<P.size(); ++i) {

    if (!B.contain(P[i])) { // process a conflict, ow keep going
      size_t ii = findPivot<dim>(P, B, i+1);
      if (ii != -1) std::swap(P[ii], P[i]);

      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();

      support.push_back(P[i]);

      // compute a ball on the prefix
      B = welzlMtfPivotSerial(P.cut(0, i), support, B);
      support.pop_back();

      if (i > dim-support.size()) {
	std::swap(P[dim-support.size()], P[i]);}
    }

  }

  return B;
}

template<int dim>
pargeo::seb::ball<dim> welzlMtfPivotParallel(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					     parlay::sequence<pargeo::point<dim>>& support,
					     pargeo::seb::ball<dim> B,
					     parlay::sequence<size_t>* flag=NULL) {
  using ballT = pargeo::seb::ball<dim>;
  using pointT = pargeo::point<dim>;

  if (P.size() < 2000) return welzlMtfPivotSerial<dim>(P, support, B);

  B = support2Ball(P, support);
  if (B.size() == dim+1) return B;

  bool freeFlag = false;
  if (!flag) {
    freeFlag = true;
    flag = new parlay::sequence<size_t>(P.size()+1);
  }

  auto process = [&](pointT p) {
                   if (!B.contain(p)) return true;
                   else return false;
                 };

  auto cleanUp = [&](parlay::slice<pointT*, pointT*> A, size_t i) {
		   size_t ii = findPivot<dim>(P, B, i+1);
		   if (ii != -1) std::swap(P[ii], P[i]);

                   if (support.size() == B.size()) B.grow(A[i]);
                   else B = ballT();
                   support.push_back(A[i]);
                   B = welzlMtfPivotParallel(A.cut(0, i), support, B, flag);
                   support.pop_back();

		   if (i > dim-support.size()) {
		     std::swap(P[dim-support.size()], P[i]);}
                 };

  parallel_prefix(P, process, cleanUp, flag, 200000, 500000);

  if(freeFlag) delete flag;
  return B;

}
