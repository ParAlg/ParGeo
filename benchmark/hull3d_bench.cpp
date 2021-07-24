#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/pseudo/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/divideConquer/hull.h"
#include "convexHull3d/parallelQuickHull/hull.h"

#include <benchmark/benchmark.h>

#include "dataLoader.h"

long maxThreads = 72;

static void serialQuick_inSphere_3d_10m(benchmark::State& state) {
  auto P = inSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_onSphere_3d_10m(benchmark::State& state) {
  auto P = onSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_inCube_3d_10m(benchmark::State& state) {
  auto P = inCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_onCube_3d_10m(benchmark::State& state) {
  auto P = onCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_thaiStatue_3d_3_5m(benchmark::State& state) {
  auto P = thaiStatue_3d_3_5m_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

static void serialQuick_lucy_3d_3_14m(benchmark::State& state) {
  auto P = lucy_3d_3_14m_float();
  for (auto _ : state) pargeo::hull3d::serialQuickHull::compute(parlay::make_slice(P));}

BENCHMARK(serialQuick_inSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_onSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_inCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_onCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_lucy_3d_3_14m)->Unit(benchmark::kMillisecond);



static void quickHull_inSphere_3d_10m(benchmark::State& state) {
  auto P = inSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_onSphere_3d_10m(benchmark::State& state) {
  auto P = onSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_inCube_3d_10m(benchmark::State& state) {
  auto P = inCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_onCube_3d_10m(benchmark::State& state) {
  auto P = onCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_thaiStatue_3d_3_5m(benchmark::State& state) {
  auto P = thaiStatue_3d_3_5m_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

static void quickHull_lucy_3d_3_14m(benchmark::State& state) {
  auto P = lucy_3d_3_14m_float();
  for (auto _ : state) pargeo::hull3d::parallelQuickHull::compute(parlay::make_slice(P), state.range(0));}

BENCHMARK(quickHull_inSphere_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_onSphere_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_inCube_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_onCube_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_buddha_3d_544k)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(quickHull_lucy_3d_3_14m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);



static void divideConquer_inSphere_3d_10m(benchmark::State& state) {
  auto P = inSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_onSphere_3d_10m(benchmark::State& state) {
  auto P = onSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_inCube_3d_10m(benchmark::State& state) {
  auto P = inCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_onCube_3d_10m(benchmark::State& state) {
  auto P = onCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_thaiStatue_3d_3_5m(benchmark::State& state) {
  auto P = thaiStatue_3d_3_5m_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

static void divideConquer_lucy_3d_3_14m(benchmark::State& state) {
  auto P = lucy_3d_3_14m_float();
  for (auto _ : state) pargeo::hull3d::divideConquer::compute(parlay::make_slice(P), state.range(0));}

BENCHMARK(divideConquer_inSphere_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_onSphere_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_inCube_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_onCube_3d_10m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_buddha_3d_544k)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);
BENCHMARK(divideConquer_lucy_3d_3_14m)->Unit(benchmark::kMillisecond)->Arg(maxThreads);



static void pseudo_inSphere_3d_10m(benchmark::State& state) {
  auto P = inSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_onSphere_3d_10m(benchmark::State& state) {
  auto P = onSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_inCube_3d_10m(benchmark::State& state) {
  auto P = inCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_onCube_3d_10m(benchmark::State& state) {
  auto P = onCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_thaiStatue_3d_3_5m(benchmark::State& state) {
  auto P = thaiStatue_3d_3_5m_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

static void pseudo_lucy_3d_3_14m(benchmark::State& state) {
  auto P = lucy_3d_3_14m_float();
  for (auto _ : state) pargeo::hull3d::pseudo::compute(parlay::make_slice(P));}

BENCHMARK(pseudo_inSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_onSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_inCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_onCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_lucy_3d_3_14m)->Unit(benchmark::kMillisecond);



static void sampling_inSphere_3d_10m(benchmark::State& state) {
  auto P = inSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_onSphere_3d_10m(benchmark::State& state) {
  auto P = onSphere_float<3>(N);
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_inCube_3d_10m(benchmark::State& state) {
  auto P = inCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_onCube_3d_10m(benchmark::State& state) {
  auto P = onCube_float<3>(N);
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_thaiStatue_3d_3_5m(benchmark::State& state) {
  auto P = thaiStatue_3d_3_5m_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

static void sampling_lucy_3d_3_14m(benchmark::State& state) {
  auto P = lucy_3d_3_14m_float();
  for (auto _ : state) pargeo::hull3d::sampling::compute(parlay::make_slice(P));}

BENCHMARK(sampling_inSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_onSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_inCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_onCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_lucy_3d_3_14m)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
