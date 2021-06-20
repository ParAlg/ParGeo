#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "pargeo/kdTree.h"
#include "pargeo/kdTreeKnn.h"
#include "pargeo/kdTreeRange.h"

static parlay::sequence<pargeo::point<2>> data2d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<2>(n, 1);
  return std::move(P);
}

static parlay::sequence<pargeo::point<3>> data3d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<3>(n, 1);
  return std::move(P);
}

static void BM_knn_2d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto P = data2d(state.range(0));
  for (auto _ : state) {
    kdNode<2, point<2>>* tree =
      buildKdt<2, point<2>>(P, true, true);
    auto nn = kdTreeKnn(P, state.range(1), tree, false);
  }
}

static void BM_range_2d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto P = data2d(state.range(0));
  for (auto _ : state) {
    kdNode<2, point<2>>* tree =
      buildKdt<2, point<2>>(P, true, true);
    parallel_for(0, P.size(), [&](size_t i) {
				auto I = kdTreeRange(P, tree, P[i], state.range(1));
			      });
  }
}

static void BM_knn_3d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto P = data3d(state.range(0));
  for (auto _ : state) {
    kdNode<3, point<3>>* tree =
      buildKdt<3, point<3>>(P, true, true);
    auto nn = kdTreeKnn(P, state.range(1), tree, false);
  }
}

static void BM_range_3d_includeTree(benchmark::State& state) {
  using namespace pargeo;
  auto P = data3d(state.range(0));
  for (auto _ : state) {
    kdNode<3, point<3>>* tree =
      buildKdt<3, point<3>>(P, true, true);
    parallel_for(0, P.size(), [&](size_t i) {
				auto I = kdTreeRange(P, tree, P[i], state.range(1));
			      });
  }
}

BENCHMARK(BM_knn_2d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Args({100000, 1})->Args({100000, 2})->Args({100000, 10})
->Args({1000000, 1})->Args({1000000, 2})->Args({1000000, 10});

BENCHMARK(BM_range_2d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000);

BENCHMARK(BM_knn_3d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Args({100000, 1})->Args({100000, 2})->Args({100000, 10})
->Args({1000000, 1})->Args({1000000, 2})->Args({1000000, 10});

BENCHMARK(BM_range_3d_includeTree)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000);

BENCHMARK_MAIN();
