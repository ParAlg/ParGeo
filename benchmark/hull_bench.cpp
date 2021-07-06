#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/pseudo/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/divideConquer/hull.h"
#include "convexHull3d/parallelQuickHull/hull.h"

#include <benchmark/benchmark.h>

#include "dataset/uniform.h"

long maxThreads = 12;
long defaultN = 1000000;

// in sphere
parlay::sequence<pargeo::fpoint<3>> data0(size_t n) {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(n, 0, 100);
}

// on sphere
parlay::sequence<pargeo::fpoint<3>> data1(size_t n) {
  return pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(n, 0, 0.05, 100);
}

// in cube
parlay::sequence<pargeo::fpoint<3>> data2(size_t n) {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(n, 1, 100);
}

static void serial_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));
}

static void serial_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));
}

static void serial_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));
}

static void parallel_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(1));
}

static void parallel_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(1));
}

static void parallel_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(1));
}

static void divide_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(1));
}

static void divide_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(1));
}

static void divide_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(1));
}

static void pseudo_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::hull3d::pseudo::compute(make_slice(P));
}

static void pseudo_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::hull3d::pseudo::compute(make_slice(P));
}

static void pseudo_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::hull3d::pseudo::compute(make_slice(P));
}

static void sampling_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::hull3d::sampling::compute(make_slice(P));
}

static void sampling_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::hull3d::sampling::compute(make_slice(P));
}

static void sampling_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::hull3d::sampling::compute(make_slice(P));
}

BENCHMARK(serial_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(serial_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(serial_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(parallel_inSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(parallel_onSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(parallel_inCube)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});

BENCHMARK(divide_inSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(divide_onSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(divide_inCube)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});

BENCHMARK(pseudo_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(pseudo_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(pseudo_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(sampling_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK_MAIN();
