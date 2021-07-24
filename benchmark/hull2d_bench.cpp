#include "convexHull2d/quickHull/hull.h"
#include "convexHull2d/divideConquer/hull.h"

#include <benchmark/benchmark.h>

#include "dataLoader.h"

static void quickHull_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere_float<2>(N);
  for (auto _ : state) pargeo::hull2d::quickHull::compute(parlay::make_slice(P));}

static void quickHull_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere_float<2>(N);
  for (auto _ : state) pargeo::hull2d::quickHull::compute(parlay::make_slice(P));}

static void quickHull_inCube_2d_10m(benchmark::State& state) {
  auto P = inCube_float<2>(N);
  for (auto _ : state) pargeo::hull2d::quickHull::compute(parlay::make_slice(P));}

static void quickHull_onCube_2d_10m(benchmark::State& state) {
  auto P = onCube_float<2>(N);
  for (auto _ : state) pargeo::hull2d::quickHull::compute(parlay::make_slice(P));}

BENCHMARK(quickHull_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(quickHull_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(quickHull_inCube_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(quickHull_onCube_2d_10m)->Unit(benchmark::kMillisecond);



static void divideConquer_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere_float<2>(N);
  for (auto _ : state) pargeo::hull2d::divideConquer::compute(parlay::make_slice(P));}

static void divideConquer_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere_float<2>(N);
  for (auto _ : state) pargeo::hull2d::divideConquer::compute(parlay::make_slice(P));}

static void divideConquer_inCube_2d_10m(benchmark::State& state) {
  auto P = inCube_float<2>(N);
  for (auto _ : state) pargeo::hull2d::divideConquer::compute(parlay::make_slice(P));}

static void divideConquer_onCube_2d_10m(benchmark::State& state) {
  auto P = onCube_float<2>(N);
  for (auto _ : state) pargeo::hull2d::divideConquer::compute(parlay::make_slice(P));}

BENCHMARK(divideConquer_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(divideConquer_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(divideConquer_inCube_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(divideConquer_onCube_2d_10m)->Unit(benchmark::kMillisecond);


BENCHMARK_MAIN();
