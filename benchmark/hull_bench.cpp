#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "convexHull3d/hull.h"

static parlay::sequence<pargeo::fpoint<3>> P[12];

/* types:
   - 0: sphere-in
   - 1: sphere-on thickness 0.05
   - 2: cube-in
   sizes: 100000/1m/10m/100m
 */
parlay::sequence<pargeo::fpoint<3>> data(size_t type, size_t size) {
  size_t i = type * 4 + size;
  if (type == 0 && size == 0) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(100000, 0));
  }
  else if (type == 0 && size == 1) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(1000000, 0));
  }
  else if (type == 0 && size == 2) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(10000000, 0));
  }
  else if (type == 0 && size == 3) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(100000000, 0));
  }
  else if (type == 1 && size == 0) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(100000, 0, 0.05));
  }
  else if (type == 1 && size == 1) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(1000000, 0, 0.05));
  }
  else if (type == 1 && size == 2) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(10000000, 0, 0.05));
  }
  else if (type == 1 && size == 3) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(100000000, 0, 0.05));
  }
  else if (type == 2 && size == 0) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(100000, 1));
  }
  else if (type == 2 && size == 1) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(1000000, 1));
  }
  else if (type == 2 && size == 2) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(10000000, 1));
  }
  else if (type == 2 && size == 3) {
    if (P[i].size() == 0)
      P[i] = std::move(pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(100000000, 1));
  }
  else {
    throw std::runtime_error("test data not implemented yet");
  }
  return P[i];
}

static void BM_hull3dSerial_inSphere_100k(benchmark::State& state) {
  auto P = data(0, 0);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_onSphere_100k(benchmark::State& state) {
  auto P = data(1, 0);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_inCube_100k(benchmark::State& state) {
  auto P = data(2, 0);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_inSphere_1m(benchmark::State& state) {
  auto P = data(0, 1);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_onSphere_1m(benchmark::State& state) {
  auto P = data(1, 1);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_inCube_1m(benchmark::State& state) {
  auto P = data(2, 1);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_inSphere_10m(benchmark::State& state) {
  auto P = data(0, 2);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_onSphere_10m(benchmark::State& state) {
  auto P = data(1, 2);
  for (auto _ : state)
    hull3dSerial(P);
}

static void BM_hull3dSerial_inCube_10m(benchmark::State& state) {
  auto P = data(2, 2);
  for (auto _ : state)
    hull3dSerial(P);
}

BENCHMARK(BM_hull3dSerial_inSphere_100k)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_onSphere_100k)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inCube_100k)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_onSphere_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inCube_1m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inSphere_10m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_onSphere_10m)->UseRealTime()->Unit(benchmark::kMillisecond);
BENCHMARK(BM_hull3dSerial_inCube_10m)->UseRealTime()->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
