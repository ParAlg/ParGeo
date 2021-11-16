#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/pseudo/hull.h"

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

BENCHMARK(serialQuick_inSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_onSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_inCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_onCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond);
BENCHMARK(serialQuick_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK(pseudo_inSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_onSphere_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_inCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_onCube_3d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_thaiStatue_3d_3_5m)->Unit(benchmark::kMillisecond);
BENCHMARK(pseudo_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
