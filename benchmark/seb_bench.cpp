#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/scan/seb.h"
#include "enclosingBall/sampling/seb.h"

#include <benchmark/benchmark.h>

#include "dataset/uniform.h"

long defaultN = 1000000;

// in sphere
parlay::sequence<pargeo::point<3>> data0(size_t n) {
  return pargeo::uniformInPolyPoints<3, pargeo::point<3>>(n, 0, 100);
}

// on sphere
parlay::sequence<pargeo::point<3>> data1(size_t n) {
  return pargeo::uniformOnPolyPoints<3, pargeo::point<3>>(n, 0, 0.05, 100);
}

// in cube
parlay::sequence<pargeo::point<3>> data2(size_t n) {
  return pargeo::uniformInPolyPoints<3, pargeo::point<3>>(n, 1, 100);
}

static void welzl_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));
}

static void welzl_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));
}

static void welzl_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));
}

static void welzlMtf_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));
}

static void welzlMtf_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));
}

static void welzlMtf_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));
}

static void welzlMtfPivot_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));
}

static void welzlMtfPivot_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));
}

static void welzlMtfPivot_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));
}

static void scan_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));
}

static void scan_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));
}

static void scan_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));
}

static void sampling_inSphere(benchmark::State& state) {
  auto P = data0(state.range(0));
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));
}

static void sampling_onSphere(benchmark::State& state) {
  auto P = data1(state.range(0));
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));
}

static void sampling_inCube(benchmark::State& state) {
  auto P = data2(state.range(0));
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));
}

BENCHMARK(welzl_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzl_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzl_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(welzlMtf_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzlMtf_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzlMtf_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(welzlMtfPivot_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzlMtfPivot_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(welzlMtfPivot_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(scan_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(scan_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(scan_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK(sampling_inSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_onSphere)->Unit(benchmark::kMillisecond)->Arg(defaultN);
BENCHMARK(sampling_inCube)->Unit(benchmark::kMillisecond)->Arg(defaultN);

BENCHMARK_MAIN();
