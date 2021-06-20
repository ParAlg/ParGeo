/* This data generator is partially adapted from the
   Problem Based Benchmark Suite
   https://github.com/cmuparlay/pbbsbench
*/

#pragma once

#include <limits>
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "pargeo/point.h"

namespace pargeo {
  namespace uniformDataGen {

    template<typename floatT>
    floatT randFloat(size_t i) {
        return floatT(parlay::hash64(i)) / floatT(std::numeric_limits<size_t>::max());
    }

    template<int dim, class pointT>
    pointT randNd(size_t i) {
			       size_t s[dim];
			       s[0] = i;
			       for (int j=1; j<dim; ++j) {
				 s[j] = j*i + parlay::hash64(s[j-1]);
			       }
			       typename pointT::floatT ss[dim];
			       for (int j=0; j<dim; ++j) {
				 ss[j] = 2 * randFloat<typename pointT::floatT>(s[j]) - 1;
			       }
			       return pointT(ss);
    }

    template<int dim, class pointT>
    pointT randInUnitSphere(size_t i) {
      auto origin = pointT();
      for(int j=0; j<dim; ++j) origin[j] = 0;
      size_t j = 0;
      pointT p;
      do {
	size_t o = parlay::hash64(j++);
	p = randNd<dim, pointT>(o+i);
      } while (p.dist(origin) > 1.0);
      return p;
    }

    template<int dim, class pointT>
    pointT randOnUnitSphere(size_t i,
				typename pointT::floatT scale=1) {
      auto origin = pointT();
      for(int j=0; j<dim; ++j) origin[j] = 0;
      pointT v = randInUnitSphere<dim, pointT>(i);
      return (v / v.dist(origin)) * scale;
    }

  } // End uniform data gen namespace

  template<int dim, class pointT = point<dim>>
  parlay::sequence<pointT> uniformInPolyPoints(size_t n,
					       size_t shape,
					       double scale = 1.0) {
    using namespace parlay;
    using namespace uniformDataGen;

    auto P = sequence<pointT>(n);
    parallel_for (0, n, [&](size_t i) {
	if (shape == 0) P[i] = randInUnitSphere<dim, pointT>(i) * scale;
	else if (shape == 1) P[i] = randNd<dim, pointT>(i) * scale;
	else throw std::runtime_error("generator not implemented yet");
      });

    return P; // data should be already permuted
  }

  template<int dim, class pointT = point<dim>>
  parlay::sequence<pointT> uniformOnPolyPoints(size_t n,
					       size_t shape,
					       double thickness,
					       double scale = 1.0) {
    using namespace parlay;
    using namespace uniformDataGen;
    using floatT = typename pointT::floatT;

    auto P = sequence<pointT>(n);

    if (shape == 0) {
      floatT r1 = 1 + thickness;
      floatT r2 = 1 - thickness;
      floatT a1 = 1; for (int d = 0; d < dim - 1; ++ d) a1 *= r1;
      floatT a2 = 1; for (int d = 0; d < dim - 1; ++ d) a2 *= r2;
      size_t n1 = a1 * n / (a1 + a2);
      size_t n2 = n - n1;
      floatT t1 = 1 - 1 / r1;
      floatT t2 = 1 / r2 - 1;

      // Outer
      parallel_for (0, n1, [&](size_t i) {
          floatT s = 1 - t1 * randFloat<floatT>(i);
	  P[i] = randOnUnitSphere<dim, pointT>(i, r1) * s * scale;
	});

      // Inner
      parallel_for (n1, n, [&](size_t i) {
          floatT s = t2 * randFloat<floatT>(i) + 1;
	  P[i] = randOnUnitSphere<dim, pointT>(i, r2) * s * scale;
	});

    } else throw std::runtime_error("generator not implemented yet");

    return P; // data should be already permuted
  }

} // End namespace
