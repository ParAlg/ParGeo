#pragma once

#include <limits>
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "pargeo/point.h"

namespace pargeo::seedSpreader {

  namespace internal {

    double randomDouble(double lo, double hi);

    template <int dim>
    class spreader : public pargeo::point<dim> {
      protected:
      using baseT = pargeo::point<dim>;

      static constexpr double lo = 0.0;
      static constexpr double hi = 100000.0;

      size_t generated;

      public:

      double rShift;
      double rVincinity;

      int cReset;
      double rhoRestart;
      int counter;

      spreader(int _cReset, double _rhoRestart);
      spreader(int _cReset, double _rhoRestart, double _rShift, double _rVincinity);

      void restart();

      pargeo::point<dim> _next();

      virtual pargeo::point<dim> next() = 0;

      void shift();

      pargeo::point<dim> random();
    };

    template <int dim>
    class simdenSpreader : public spreader<dim> {
      using baseT = spreader<dim>;
      public:
      pargeo::point<dim> next();
      simdenSpreader(int _cReset, double _rhoRestart, double _rShift, double rVincinity);
    };

    template <int dim>
    class vardenSpreader : public spreader<dim> {
      private:
      using baseT = spreader<dim>;
      public:
      void step();
      pargeo::point<dim> next();
      vardenSpreader(int cReset, double rhoRestart);
    };

  } // End namespace pargeo::seedSpreader::internal

  template <int dim, class pointT = pargeo::point<dim>>
  parlay::sequence<pointT> simdenGenerator(size_t n, double rhoNoise = 0.0001);

  template <int dim, class pointT = pargeo::point<dim>>
  parlay::sequence<pointT> vardenGenerator(size_t n, double rhoNoise = 0.0001);

} // End namespace pargeo::seedSpreader
