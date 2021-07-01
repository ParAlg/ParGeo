#include "convexHull3d/serialHull.h"
#include "convexHull3d/pseudoHull.h"
#include "convexHull3d/gridHull.h"
#include "convexHull3d/samplingHull.h"
#include "convexHull3d/searchHull.h"
#include "convexHull3d/concurrentHull.h"
// #include "convexHull3d/incrementalHull.h"
#include "convexHull3d/parallelHull.h"
#include "convexHull3d/samplingHull.h"
#include "convexHull3d/searchHull.h"

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
  for (auto _ : state) hull3dSerial(P);
}

static void serial_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dSerial(P);
}

static void serial_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dSerial(P);
}

// static void incremental_inSphere(benchmark::State& state) {
//   auto P = data0(state.range(0));
//   for (auto _ : state) hull3dIncremental(P, state.range(1));
// }

// static void incremental_onSphere(benchmark::State& state) {
//   auto P = data1(state.range(0));
//   for (auto _ : state) hull3dIncremental(P, state.range(1));
// }

// static void incremental_inCube(benchmark::State& state) {
//   auto P = data2(state.range(0));
//   for (auto _ : state) hull3dIncremental(P, state.range(1));
// }

static void parallel_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dParallel(P, state.range(1));
}

static void parallel_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dParallel(P, state.range(1));
}

static void parallel_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dParallel(P, state.range(1));
}

static void concurrent_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dConcurrent(P, state.range(1));
}

static void concurrent_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dConcurrent(P, state.range(1));
}

static void concurrent_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dConcurrent(P, state.range(1));
}

static void grid_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void grid_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void grid_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void pseudo_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dPseudo(P);
}

static void pseudo_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dPseudo(P);
}

static void pseudo_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dPseudo(P);
}

static void sampling_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dSampling(P);
}

static void sampling_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dSampling(P);
}

static void sampling_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dSampling(P);
}

static void search_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) hull3dSearch(P);
}

static void search_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) hull3dSearch(P);
}

static void search_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) hull3dSearch(P);
}

BENCHMARK(serial_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(serial_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(serial_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

// BENCHMARK(incremental_inSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
// BENCHMARK(incremental_onSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
// BENCHMARK(incremental_inCube)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});

BENCHMARK(parallel_inSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(parallel_onSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(parallel_inCube)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});

BENCHMARK(concurrent_inSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(concurrent_onSphere)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});
BENCHMARK(concurrent_inCube)->Unit(benchmark::kMillisecond)->Args({defaultN, maxThreads});

BENCHMARK(grid_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(grid_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(grid_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(pseudo_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(pseudo_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(pseudo_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(sampling_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(search_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(search_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(search_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK_MAIN();
