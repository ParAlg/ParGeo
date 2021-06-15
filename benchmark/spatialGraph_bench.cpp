#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "spatialGraph/spatialGraph.h"

static parlay::sequence<pargeo::point<2>> P;

void data() {
  if (P.size() == 0)
    P = pargeo::uniformInPolyPoints<2>(1000000, 1);
}

static void BM_knnGraph(benchmark::State& state) {
  data();
  for (auto _ : state)
    parlay::sequence<pargeo::dirEdge> E = pargeo::knnGraph(P, state.range(0));
}

static void BM_delaunayGraph(benchmark::State& state) {
  data();
  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::delaunayGraph(P);
}

static void BM_gabrielGraph(benchmark::State& state) {
  data();
  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::gabrielGraph(P);
}

static void BM_betaSkeleton(benchmark::State& state) {
  data();
  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::betaSkeleton(P, state.range(0));
}

static void BM_spanner(benchmark::State& state) {
  data();
  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::spanner(P, state.range(0));
}

BENCHMARK(BM_knnGraph)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(1)->Arg(3)->Arg(10);

BENCHMARK(BM_delaunayGraph)
->UseRealTime()
->Unit(benchmark::kMillisecond);

BENCHMARK(BM_gabrielGraph)
->UseRealTime()
->Unit(benchmark::kMillisecond);

BENCHMARK(BM_betaSkeleton)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(1)->Arg(3)->Arg(10);

BENCHMARK(BM_spanner)
->UseRealTime()
->Unit(benchmark::kMillisecond)
->Arg(10)->Arg(20)->Arg(30);

BENCHMARK_MAIN();
