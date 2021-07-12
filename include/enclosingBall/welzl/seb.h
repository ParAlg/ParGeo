#pragma once

#include "parlay/sequence.h"
#include "pargeo/point.h"
#include "enclosingBall/ball.h"

namespace pargeo {
  namespace seb {

    namespace welzl {
      template<int dim>
      pargeo::seb::ball<dim> compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*>);
    }

    namespace welzlMtf {
      template<int dim>
      pargeo::seb::ball<dim> compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*>);
    }

    namespace welzlMtfPivot {
      template<int dim>
      pargeo::seb::ball<dim> compute(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*>);
    }

  }
}
