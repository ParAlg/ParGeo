#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "convexHull3d/hull.h"

static parlay::sequence<pargeo::fpoint<3>> P[3];

/* types:
   - 0: sphere-in
   - 1: sphere-on thickness 0.05
   - 2: cube-in
 */
parlay::sequence<pargeo::fpoint<3>> data(size_t type) {
  size_t i = type;
  if (type == 0) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(1000000, 0));
  }
  else if (type == 1) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(1000000, 0, 0.05));
  }
  else if (type == 2) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(1000000, 1));
  }
  else {
    throw std::runtime_error("test data not implemented yet");
  }
  return P[i];
}

static void BM_hull3dSerial_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dSerial(P);
}

static void BM_hull3dSerial_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dSerial(P);
}

static void BM_hull3dSerial_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dSerial(P);
}

static void BM_hull3dIncremental_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dIncremental(P, state.range(0));
}

static void BM_hull3dIncremental_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dIncremental(P, state.range(0));
}

static void BM_hull3dIncremental_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dIncremental(P, state.range(0));
}

static void BM_hull3dConcurrent_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dConcurrent(P, state.range(0));
}

static void BM_hull3dConcurrent_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dConcurrent(P, state.range(0));
}

static void BM_hull3dConcurrent_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dConcurrent(P, state.range(0));
}

static void BM_hull3dGrid_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void BM_hull3dGrid_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void BM_hull3dGrid_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dGrid(P, 4, false);
}

static void BM_hull3dPseudo_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dPseudo(P);
}

static void BM_hull3dPseudo_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dPseudo(P);
}

static void BM_hull3dPseudo_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dPseudo(P);
}

static void BM_hull3dSampling_inSphere_1m(benchmark::State& state) {
  auto P = data(0);
  for (auto _ : state) hull3dSampling(P, 0.01);
}

static void BM_hull3dSampling_onSphere_1m(benchmark::State& state) {
  auto P = data(1);
  for (auto _ : state) hull3dSampling(P, 0.01);
}

static void BM_hull3dSampling_inCube_1m(benchmark::State& state) {
  auto P = data(2);
  for (auto _ : state) hull3dSampling(P, 0.01);
}

BENCHMARK(BM_hull3dSerial_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond);

BENCHMARK(BM_hull3dIncremental_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);
BENCHMARK(BM_hull3dIncremental_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);
BENCHMARK(BM_hull3dIncremental_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);

BENCHMARK(BM_hull3dConcurrent_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);
BENCHMARK(BM_hull3dConcurrent_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);
BENCHMARK(BM_hull3dConcurrent_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond)->Arg(12)->Arg(1);

BENCHMARK(BM_hull3dGrid_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dGrid_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dGrid_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond);

BENCHMARK(BM_hull3dPseudo_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dPseudo_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dPseudo_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
