#include "enclosingBall/scan/seb.h"
#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/sampling/seb.h"

#include <benchmark/benchmark.h>

#include "dataLoader.h"

static void welzl_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzl::compute<2>(parlay::make_slice(P));}

static void welzl_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzl::compute<2>(parlay::make_slice(P));}

static void welzl_inSphere_5d_10m(benchmark::State& state) {
  auto P = inSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzl::compute<5>(parlay::make_slice(P));}

static void welzl_onSphere_5d_10m(benchmark::State& state) {
  auto P = onSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzl::compute<5>(parlay::make_slice(P));}

static void welzl_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k();
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));}

static void welzl_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k();
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));}

static void welzl_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k();
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));}

static void welzl_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m();
  for (auto _ : state) pargeo::seb::welzl::compute<3>(parlay::make_slice(P));}



static void welzlMtf_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzlMtf::compute<2>(parlay::make_slice(P));}

static void welzlMtf_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzlMtf::compute<2>(parlay::make_slice(P));}

static void welzlMtf_inSphere_5d_10m(benchmark::State& state) {
  auto P = inSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzlMtf::compute<5>(parlay::make_slice(P));}

static void welzlMtf_onSphere_5d_10m(benchmark::State& state) {
  auto P = onSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzlMtf::compute<5>(parlay::make_slice(P));}

static void welzlMtf_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k();
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));}

static void welzlMtf_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k();
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));}

static void welzlMtf_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k();
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));}

static void welzlMtf_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m();
  for (auto _ : state) pargeo::seb::welzlMtf::compute<3>(parlay::make_slice(P));}



static void welzlMtfPivot_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<2>(parlay::make_slice(P));}

static void welzlMtfPivot_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere<2>(N);
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<2>(parlay::make_slice(P));}

static void welzlMtfPivot_inSphere_5d_10m(benchmark::State& state) {
  auto P = inSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<5>(parlay::make_slice(P));}

static void welzlMtfPivot_onSphere_5d_10m(benchmark::State& state) {
  auto P = onSphere<5>(N);
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<5>(parlay::make_slice(P));}

static void welzlMtfPivot_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k();
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));}

static void welzlMtfPivot_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k();
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));}

static void welzlMtfPivot_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k();
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));}

static void welzlMtfPivot_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m();
  for (auto _ : state) pargeo::seb::welzlMtfPivot::compute<3>(parlay::make_slice(P));}



static void scan_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere<2>(N);
  for (auto _ : state) pargeo::seb::scan::compute<2>(parlay::make_slice(P));}

static void scan_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere<2>(N);
  for (auto _ : state) pargeo::seb::scan::compute<2>(parlay::make_slice(P));}

static void scan_inSphere_5d_10m(benchmark::State& state) {
  auto P = inSphere<5>(N);
  for (auto _ : state) pargeo::seb::scan::compute<5>(parlay::make_slice(P));}

static void scan_onSphere_5d_10m(benchmark::State& state) {
  auto P = onSphere<5>(N);
  for (auto _ : state) pargeo::seb::scan::compute<5>(parlay::make_slice(P));}

static void scan_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k();
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));}

static void scan_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k();
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));}

static void scan_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k();
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));}

static void scan_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m();
  for (auto _ : state) pargeo::seb::scan::compute<3>(parlay::make_slice(P));}



static void sampling_inSphere_2d_10m(benchmark::State& state) {
  auto P = inSphere<2>(N);
  for (auto _ : state) pargeo::seb::sampling::compute<2>(parlay::make_slice(P));}

static void sampling_onSphere_2d_10m(benchmark::State& state) {
  auto P = onSphere<2>(N);
  for (auto _ : state) pargeo::seb::sampling::compute<2>(parlay::make_slice(P));}

static void sampling_inSphere_5d_10m(benchmark::State& state) {
  auto P = inSphere<5>(N);
  for (auto _ : state) pargeo::seb::sampling::compute<5>(parlay::make_slice(P));}

static void sampling_onSphere_5d_10m(benchmark::State& state) {
  auto P = onSphere<5>(N);
  for (auto _ : state) pargeo::seb::sampling::compute<5>(parlay::make_slice(P));}

static void sampling_armadillo_3d_173k(benchmark::State& state) {
  auto P = armadillo_3d_173k();
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));}

static void sampling_dragon_3d_438k(benchmark::State& state) {
  auto P = dragon_3d_438k();
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));}

static void sampling_buddha_3d_544k(benchmark::State& state) {
  auto P = buddha_3d_544k();
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));}

static void sampling_asianDragon_3d_3_6m(benchmark::State& state) {
  auto P = asianDragon_3d_3_6m();
  for (auto _ : state) pargeo::seb::sampling::compute<3>(parlay::make_slice(P));}

// BENCHMARK(welzl_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_inSphere_5d_10m)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_onSphere_5d_10m)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_armadillo_3d_173k)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_dragon_3d_438k)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_buddha_3d_544k)->Unit(benchmark::kMillisecond);
// BENCHMARK(welzl_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK(welzlMtf_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_inSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_onSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_armadillo_3d_173k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtf_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK(welzlMtfPivot_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_inSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_onSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_armadillo_3d_173k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(welzlMtfPivot_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK(scan_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_inSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_onSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_armadillo_3d_173k)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(scan_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK(sampling_inSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_onSphere_2d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_inSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_onSphere_5d_10m)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_armadillo_3d_173k)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_dragon_3d_438k)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_buddha_3d_544k)->Unit(benchmark::kMillisecond);
BENCHMARK(sampling_asianDragon_3d_3_6m)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
