#include "parlay/random.h"
#include "parlay/sequence.h"
#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/welzl/welzl.h"
#include "enclosingBall/ball.h"

//todo make point generic

template<int dim>
pargeo::seb::ball<dim>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> _P) {

  typedef pargeo::point<dim> pointT;
  typedef pargeo::seb::ball<dim> ballT;

  auto support = parlay::sequence<pointT>();
  auto P = parlay::random_shuffle(_P);
  ballT D = pargeo::seb::welzl::welzlMtfPivotParallel(parlay::make_slice(P), support, ballT());
  return D;
}

template
pargeo::seb::ball<2>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*> P);

template
pargeo::seb::ball<3>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*> P);

template
pargeo::seb::ball<4>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<4>*, pargeo::point<4>*> P);

template
pargeo::seb::ball<5>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<5>*, pargeo::point<5>*> P);

template
pargeo::seb::ball<6>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<6>*, pargeo::point<6>*> P);

template
pargeo::seb::ball<7>
pargeo::seb::welzlMtfPivot::compute(parlay::slice<pargeo::point<7>*, pargeo::point<7>*> P);
