/* The seed spreader is based on these papers:
[1] Junhao Gan, Yufei Tao. On the Hardness and Approximation of Euclidean DBSCAN
[2] Junhao Gan, Yufei Tao. DBSCAN Revisited: Mis-Claim, Un-Fixability, and Approximation
*/

#include "dataset/seedSpreader.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"

double pargeo::seedSpreader::internal::randomDouble(double lo, double hi)
{
  double f = (double) rand() / RAND_MAX;
  return lo + f * (hi - lo);
}

/* spreader base class implementations */

template <int dim>
pargeo::seedSpreader::internal::spreader<dim>::spreader(
  int _cReset,
  double _rhoRestart):
  pargeo::point<dim>(),
  counter(0),
  cReset(_cReset),
  rhoRestart(_rhoRestart),
  generated(0)
{
  restart();
}

template <int dim>
pargeo::seedSpreader::internal::spreader<dim>::spreader(
  int _cReset,
  double _rhoRestart,
  double _rShift,
  double _rVincinity):
  pargeo::point<dim>(),
  counter(0),
  cReset(_cReset),
  rhoRestart(_rhoRestart),
  rShift(_rShift),
  rVincinity(_rVincinity),
  generated(0)
{
  restart();
}

template <int dim>
void pargeo::seedSpreader::internal::spreader<dim>::restart() {
  counter = 0;
  for (int i = 0; i < dim; ++ i) baseT::x[i] = randomDouble(lo, hi);
}

template <int dim>
void pargeo::seedSpreader::internal::spreader<dim>::shift() {
  counter = 0;

  // generate a random direction
  pargeo::point<dim> dir = pargeo::point<dim>();
  for (int i = 0; i < dim; ++ i) dir.x[i] = randomDouble(-1, 1);

  // scale the vector to the desired length
  double len = dir.length();
  for (int i = 0; i < dim; ++ i) baseT::x[i] += rShift * dir.x[i] / len;
}

template <int dim>
pargeo::point<dim>
pargeo::seedSpreader::internal::spreader<dim>::random() {
  auto pt = pargeo::point<dim>();
  for (int i = 0; i < dim; ++ i) pt[i] = randomDouble(lo, hi);
  return pt;
}

template<int dim>
pargeo::point<dim>
pargeo::seedSpreader::internal::spreader<dim>::_next() {
  if (randomDouble(0, 1) <= rhoRestart) restart();
  if (counter >= cReset) shift();
  pargeo::point<dim> pt(baseT::x);
  do {
    for (int i = 0; i < dim; ++ i)
      pt.x[i] = baseT::x[i] + randomDouble(-rVincinity, rVincinity);
  } while (baseT::dist(pt) > rVincinity);
  counter ++; // local
  generated ++; // global
  return pt;
}

/* simdenSpreader class implementations */

template <int dim>
pargeo::seedSpreader::internal::simdenSpreader<dim>::simdenSpreader(int _cReset, double _rhoRestart,
  double _rShift, double _rVincinity):
  spreader<dim>(_cReset, _rhoRestart, _rShift, _rVincinity)
{
}

template<int dim>
pargeo::point<dim>
pargeo::seedSpreader::internal::simdenSpreader<dim>::next() {
  return baseT::_next();
}

/* vardenSpreader class implementations */

template <int dim>
pargeo::seedSpreader::internal::vardenSpreader<dim>::vardenSpreader(
  int _cReset, double _rhoRestart):
  spreader<dim>(_cReset, _rhoRestart)
{
  step(); // optional call
}

template <int dim>
void pargeo::seedSpreader::internal::vardenSpreader<dim>::step() {
  baseT::rVincinity = 100 * ((baseT::generated % 10) + 1);
  baseT::rShift = baseT::rVincinity * dim / 2;
}

template<int dim>
pargeo::point<dim>
pargeo::seedSpreader::internal::vardenSpreader<dim>::next() {
  step();
  return baseT::_next();
}

/* Generator implementations */

template <int dim, class pointT = pargeo::point<dim>>
parlay::sequence<pointT>
pargeo::seedSpreader::simdenGenerator(size_t n, double rhoNoise) {
  srand(123);

  auto sp = internal::simdenSpreader<dim>(
    100,
    10.0 / (n * (1.0 - rhoNoise)),
    50 * dim,
    100);

  parlay::sequence<pointT> data(n);

  size_t p = 0;

  // Generate clustered points

  for (; p < (1 - rhoNoise) * n; p ++) data[p] = sp.next();

  // Generate noise

  for (; p < n; p ++) data[p] = sp.random();

  return data;
}

template <int dim, class pointT = pargeo::point<dim>>
parlay::sequence<pointT>
pargeo::seedSpreader::vardenGenerator(size_t n, double rhoNoise) {
  srand(321);

  auto sp = internal::vardenSpreader<dim>(100, 10.0 / (n * (1.0 - rhoNoise)));

  parlay::sequence<pointT> data(n);

  size_t p = 0;

  // Generate clustered points

  for (; p < (1 - rhoNoise) * n; p ++) data[p] = sp.next();

  // Generate noise

  for (; p < n; p ++) data[p] = sp.random();

  return data;
}

// TODO add other dim definitions

template parlay::sequence<pargeo::point<2>>
pargeo::seedSpreader::simdenGenerator<2>(size_t, double);

template parlay::sequence<pargeo::point<2>>
pargeo::seedSpreader::vardenGenerator<2>(size_t, double);
