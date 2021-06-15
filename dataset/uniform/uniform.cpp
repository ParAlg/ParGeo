/* This data generator is partially adapted from the
   Problem Based Benchmark Suite
   https://github.com/cmuparlay/pbbsbench
*/

#include <limits>
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "pargeo/point.h"

#include "dataset/dataset.h"

using namespace parlay;

namespace pargeo {

using floatT = double;

double randDouble(size_t i) {
  return double(parlay::hash64(i)) / double(std::numeric_limits<size_t>::max());
}

template<int dim>
point<dim> randNd(size_t i) {
  size_t s[dim];
  s[0] = i;
  for (int j=1; j<dim; ++j) {
    s[j] = j*i + parlay::hash64(s[j-1]);
  }
  floatT ss[dim];
  for (int j=0; j<dim; ++j) {
    ss[j] = 2 * randDouble(s[j]) - 1;
  }
  return point<dim>(ss);
}

template<int dim>
point<dim> randInUnitSphere(size_t i) {
  auto origin = point<dim>();
  for(int j=0; j<dim; ++j) origin[j] = 0;
  size_t j = 0;
  point<dim> p;
  do {
    size_t o = parlay::hash64(j++);
    p = randNd<dim>(o+i);
  } while (p.dist(origin) > 1.0);
  return p;
}

template<int dim>
point<dim> randOnUnitSphere(size_t i, floatT scale=1) {
  auto origin = point<dim>();
  for(int j=0; j<dim; ++j) origin[j] = 0;
  point<dim> v = randInUnitSphere<dim>(i);
  return (v / v.dist(origin)) * scale;
}

template<int dim>
sequence<point<dim>> uniformInPolyPoints(size_t n, size_t shape, double scale) {

  auto P = sequence<point<dim>>(n);
  parallel_for (0, n, [&](size_t i) {
			if (shape == 0) P[i] = randInUnitSphere<dim>(i) * scale;
			else if (shape == 1) P[i] = randNd<dim>(i) * scale;
			else throw std::runtime_error("generator not implemented yet");
		      });

  return P; // data should be already permuted
}

template<int dim>
sequence<point<dim>> uniformOnPolyPoints(size_t n, size_t shape, double thickness, double scale) {

  auto P = sequence<point<dim>>(n);

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
			   floatT s = 1 - t1 * randDouble(i);
			   P[i] = randOnUnitSphere<dim>(i, r1) * s * scale;
			 });

    // Inner
    parallel_for (n1, n, [&](size_t i) {
			   floatT s = t2 * randDouble(i) + 1;
			   P[i] = randOnUnitSphere<dim>(i, r2) * s * scale;
			 });

  } else throw std::runtime_error("generator not implemented yet");

  return P; // data should be already permuted
}

template parlay::sequence<pargeo::point<2>> pargeo::uniformInPolyPoints<2>(size_t, size_t, double);
template parlay::sequence<pargeo::point<3>> pargeo::uniformInPolyPoints<3>(size_t, size_t, double);
template parlay::sequence<pargeo::point<4>> pargeo::uniformInPolyPoints<4>(size_t, size_t, double);
template parlay::sequence<pargeo::point<5>> pargeo::uniformInPolyPoints<5>(size_t, size_t, double);
template parlay::sequence<pargeo::point<6>> pargeo::uniformInPolyPoints<6>(size_t, size_t, double);
template parlay::sequence<pargeo::point<7>> pargeo::uniformInPolyPoints<7>(size_t, size_t, double);
template parlay::sequence<pargeo::point<8>> pargeo::uniformInPolyPoints<8>(size_t, size_t, double);
template parlay::sequence<pargeo::point<9>> pargeo::uniformInPolyPoints<9>(size_t, size_t, double);

template parlay::sequence<pargeo::point<2>> pargeo::uniformOnPolyPoints<2>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<3>> pargeo::uniformOnPolyPoints<3>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<4>> pargeo::uniformOnPolyPoints<4>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<5>> pargeo::uniformOnPolyPoints<5>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<6>> pargeo::uniformOnPolyPoints<6>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<7>> pargeo::uniformOnPolyPoints<7>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<8>> pargeo::uniformOnPolyPoints<8>(size_t, size_t, double, double);
template parlay::sequence<pargeo::point<9>> pargeo::uniformOnPolyPoints<9>(size_t, size_t, double, double);

} // End namespace
