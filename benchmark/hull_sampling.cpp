#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "convexHull3d/sampling/hull.h"

long maxThreads = 12;
long defaultN = 1000000;

// in sphere
parlay::sequence<pargeo::fpoint<3>> data0(size_t n = defaultN) {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(n, 0, 100);
}

// on sphere
parlay::sequence<pargeo::fpoint<3>> data1(size_t n = defaultN) {
  return pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(n, 0, 0.05, 100);
}

// in cube
parlay::sequence<pargeo::fpoint<3>> data2(size_t n = defaultN) {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(n, 1, 100);
}

static void randomSampling_inSphere(benchmark::State& state, double fraction) {
  auto P = data0();
  for (auto _ : state) pargeo::hull3d::sampling::random(make_slice(P), fraction);
}

static void randomSampling_onSphere(benchmark::State& state, double fraction) {
  auto P = data1();
  for (auto _ : state) pargeo::hull3d::sampling::random(make_slice(P), fraction);
}

static void randomSampling_inCube(benchmark::State& state, double fraction) {
  auto P = data2();
  for (auto _ : state) pargeo::hull3d::sampling::random(make_slice(P), fraction);
}

BENCHMARK_CAPTURE(randomSampling_inSphere, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(randomSampling_onSphere, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_onSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_onSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_onSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(randomSampling_inCube, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inCube, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inCube, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomSampling_inCube, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

static void gridSampling_inSphere(benchmark::State& state, double fraction) {
  auto P = data0();
  for (auto _ : state) pargeo::hull3d::sampling::grid(make_slice(P), fraction);
}

static void gridSampling_onSphere(benchmark::State& state, double fraction) {
  auto P = data1();
  for (auto _ : state) pargeo::hull3d::sampling::grid(make_slice(P), fraction);
}

static void gridSampling_inCube(benchmark::State& state, double fraction) {
  auto P = data2();
  for (auto _ : state) pargeo::hull3d::sampling::grid(make_slice(P), fraction);
}

BENCHMARK_CAPTURE(gridSampling_inSphere, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(gridSampling_onSphere, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_onSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_onSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_onSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(gridSampling_inCube, 0.1, 0.1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inCube, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inCube, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(gridSampling_inCube, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);

static void randomProjection_inSphere(benchmark::State& state, double fraction) {
  auto P = data0();
  for (auto _ : state) pargeo::hull3d::sampling::projection(make_slice(P), fraction);
}

static void randomProjection_onSphere(benchmark::State& state, double fraction) {
  auto P = data1();
  for (auto _ : state) pargeo::hull3d::sampling::projection(make_slice(P), fraction);
}

static void randomProjection_inCube(benchmark::State& state, double fraction) {
  auto P = data2();
  for (auto _ : state) pargeo::hull3d::sampling::projection(make_slice(P), fraction);
}

BENCHMARK_CAPTURE(randomProjection_inSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inSphere, 0.00001, 0.00001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(randomProjection_onSphere, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_onSphere, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_onSphere, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_onSphere, 0.00001, 0.00001)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(randomProjection_inCube, 0.01, 0.01)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inCube, 0.001, 0.001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inCube, 0.0001, 0.0001)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(randomProjection_inCube, 0.00001, 0.00001)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
