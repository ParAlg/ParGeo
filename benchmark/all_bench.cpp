#include <benchmark/benchmark.h>

#include "dataset/uniform.h"
#include "parlay/utilities.h"
#include "pargeo/kdTree.h"
#include "pargeo/kdTreeKnn.h"
#include "pargeo/kdTreeRange.h"
#include "pargeo/bccp.h"
#include "pargeo/wspd.h"

int N = 1000000;

static parlay::sequence<pargeo::point<2>> data2d() {
  auto P = pargeo::uniformInPolyPoints<2>(N, 1);
  return std::move(P);
}

static parlay::sequence<pargeo::point<3>> data3d() {
  auto P = pargeo::uniformInPolyPoints<3>(N, 1);
  return std::move(P);
}

static parlay::sequence<pargeo::point<5>> data5d() {
  auto P = pargeo::uniformInPolyPoints<5>(N, 1);
  return std::move(P);
}

static void kdTree_build_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  for (auto _ : state) {
    kdNode<2, point<2>>* tree =
      buildKdt<2, point<2>>(points, true, true);
  }
}

static void kdTree_knn_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdNode<2, point<2>>* tree =
    buildKdt<2, point<2>>(points, true, true);

  for (auto _ : state) {
    kdTreeKnn(points, 2, tree);
  }

}

static void kdTree_rangeSearch_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdNode<2, point<2>>* tree =
    buildKdt<2, point<2>>(points, true, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Spherical range search with random radius
      kdTreeRange(points, tree, points[i],
		  parlay::hash64(i)/numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_orthRangeSearch_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdNode<2, point<2>>* tree =
    buildKdt<2, point<2>>(points, true, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Orthogonal range search with random box
      kdTreeOrthRange(points, tree, points[i],
		      parlay::hash64(i)/numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_bccp_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdNode<2, point<2>>* tree1 =
    buildKdt<2, point<2>>(points.cut(0, N/2), true, true);

  kdNode<2, point<2>>* tree2 =
    buildKdt<2, point<2>>(points.cut(N/2, N), true, true);

  for (auto _ : state) {
    bccp(tree1, tree2);
  }
}

static void kdTree_wspd_s2_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdNode<2, point<2>>* tree =
    buildKdt<2, point<2>>(points, true, true);

  for (auto _ : state) {
    wspdParallel(tree, 2);
  }
}

BENCHMARK(kdTree_build_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_knn_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_rangeSearch_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_orthRangeSearch_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_bccp_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_wspd_s2_2d)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
