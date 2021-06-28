#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "pargeo/kdTree.h"

static parlay::sequence<pargeo::point<2>> data2d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<2>(n, 1);
  return std::move(P);
}

static parlay::sequence<pargeo::point<3>> data3d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<3>(n, 1);
  return std::move(P);
}

static void BM_kdTree_2d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d(state.range(0));
  for (auto _ : state) {
    kdNode<2, point<2>>* tree =
      buildKdt<2, point<2>>(points, true, true);
  }
}

static void BM_kdTree_3d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto points = data3d(state.range(0));
  for (auto _ : state) {
    kdNode<3, point<3>>* tree =
      buildKdt<3, point<3>>(points, true, true);
  }
}

BENCHMARK(BM_kdTree_2d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000);

BENCHMARK(BM_kdTree_3d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000);

BENCHMARK_MAIN();
