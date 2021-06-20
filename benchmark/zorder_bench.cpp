#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "pargeo/zorderSort.h"

static parlay::sequence<pargeo::point<2>> data2d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<2>(n, 1);
  return std::move(P);
}

static parlay::sequence<pargeo::point<3>> data3d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<3>(n, 1);
  return std::move(P);
}

static void BM_zorder_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d(state.range(0));
  for (auto _ : state) {
    zorderSort2d(points);
  }
}

static void BM_zorder_3d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data3d(state.range(0));
  for (auto _ : state) {
    zorderSort3d(points);
  }
}

BENCHMARK(BM_zorder_2d)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000)->Arg(10000000);

BENCHMARK(BM_zorder_3d)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(100000)->Arg(1000000)->Arg(10000000);

BENCHMARK_MAIN();
