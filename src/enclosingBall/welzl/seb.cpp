#include <vector>
#include "parlay/sequence.h"
#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/welzl/welzl.h"
#include "enclosingBall/ball.h"

//todo make point generic

template<int dim>
pargeo::seb::ball<dim>
pargeo::seb::welzl::compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P) {

  typedef pargeo::point<dim> pointT;
  typedef pargeo::seb::ball<dim> ballT;

  auto support = parlay::sequence<pointT>();
  return pargeo::seb::welzl::welzlParallel(P, support, ballT());
}

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