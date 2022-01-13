#include <benchmark/benchmark.h>

#include "closestPair/closestPair.h"
#include "convexHull2d/divideConquer/hull.h"
#include "convexHull3d/pseudo/hull.h"
#include "dataset/uniform.h"
#include "euclideanMst/euclideanMst.h"
#include "enclosingBall/sampling/seb.h"
#include "kdTree/kdTree.h"
#include "mortonSort/mortonSort.h"
#include "parlay/utilities.h"
#include "pargeo/dynKdTree.h"
#include "spatialGraph/spatialGraph.h"

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
    kdTree::node<2, point<2>>* tree =
      kdTree::build<2, point<2>>(points, true);
  }
}

static void kdTree_knn_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdTree::node<2, point<2>>* tree =
    kdTree::build<2, point<2>>(points, true);

  for (auto _ : state) {
    kdTree::batchKnn(points, 2, tree);
  }

}

static void kdTree_rangeSearch_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdTree::node<2, point<2>>* tree =
    kdTree::build<2, point<2>>(points, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Spherical range search with random radius
      kdTree::rangeSearch(points, tree, points[i],
		  parlay::hash64(i)/std::numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_orthRangeSearch_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdTree::node<2, point<2>>* tree =
    kdTree::build<2, point<2>>(points, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Orthogonal range search with random box
      kdTree::orthogonalRangeSearch(points, tree, points[i],
		      parlay::hash64(i)/std::numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_bccp_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdTree::node<2, point<2>>* tree1 =
    kdTree::build<2, point<2>>(points.cut(0, N/2), true);

  kdTree::node<2, point<2>>* tree2 =
    kdTree::build<2, point<2>>(points.cut(N/2, N), true);

  for (auto _ : state) {
    kdTree::bccp(tree1, tree2);
  }
}

static void kdTree_wspd_s2_2d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data2d();

  kdTree::node<2, point<2>>* tree =
    kdTree::build<2, point<2>>(points, true, 1);

  for (auto _ : state) {
    kdTree::wellSeparatedPairDecomp(tree, 2);
  }
}

static void kdTree_build_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  for (auto _ : state) {
    kdTree::node<5, point<5>>* tree =
      kdTree::build<5, point<5>>(points, true);
  }
}

static void kdTree_knn_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  kdTree::node<5, point<5>>* tree =
    kdTree::build<5, point<5>>(points, true);

  for (auto _ : state) {
    kdTree::batchKnn(points, 5, tree);
  }

}

static void kdTree_rangeSearch_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  kdTree::node<5, point<5>>* tree =
    kdTree::build<5, point<5>>(points, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Spherical range search with random radius
      kdTree::rangeSearch(points, tree, points[i],
		  parlay::hash64(i)/std::numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_orthRangeSearch_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  kdTree::node<5, point<5>>* tree =
    kdTree::build<5, point<5>>(points, true);

  for (auto _ : state) {

    parlay::parallel_for(0, N, [&](size_t i){

      // Orthogonal range search with random box
      kdTree::orthogonalRangeSearch(points, tree, points[i],
		      parlay::hash64(i)/std::numeric_limits<size_t>::max());

    });

  }
}

static void kdTree_bccp_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  kdTree::node<5, point<5>>* tree1 =
    kdTree::build<5, point<5>>(points.cut(0, N/2), true);

  kdTree::node<5, point<5>>* tree2 =
    kdTree::build<5, point<5>>(points.cut(N/2, N), true);

  for (auto _ : state) {
    kdTree::bccp(tree1, tree2);
  }
}

static void kdTree_wspd_s2_5d(benchmark::State& state) {
  using namespace pargeo;
  auto points = data5d();

  kdTree::node<5, point<5>>* tree =
    kdTree::build<5, point<5>>(points, true);

  for (auto _ : state) {
    kdTree::wellSeparatedPairDecomp(tree, 2);
  }
}

static void dynamic_kdTree_build_2d(benchmark::State& state) {
  using namespace pargeo;
  using namespace pargeo::dynKdTree;
  static const int dim = 2;
  auto P = data2d();

  using nodeT = rootNode<dim, pargeo::point<dim>>;

  for (auto _ : state) {
    std::unique_ptr<nodeT> tree1 = std::unique_ptr<nodeT>(new nodeT(P));
  }
}

static void dynamic_kdTree_insert_2d(benchmark::State& state) {
  using namespace pargeo;
  using namespace pargeo::dynKdTree;
  static const int dim = 2;
  auto P = data2d();

  using nodeT = rootNode<dim, pargeo::point<dim>>;

  std::unique_ptr<nodeT> tree1 = std::unique_ptr<nodeT>(new nodeT(P));

  for (auto _ : state) {
    tree1->insert(P);
  }
}

static void dynamic_kdTree_erase_2d(benchmark::State& state) {
  using namespace pargeo;
  using namespace pargeo::dynKdTree;
  static const int dim = 2;
  auto P = data2d();

  using nodeT = rootNode<dim, pargeo::point<dim>>;

  std::unique_ptr<nodeT> tree1 = std::unique_ptr<nodeT>(new nodeT(P));

  for (auto _ : state) {
    tree1->erase(P);
  }
}

static void emst_2d(benchmark::State& state) {
  using namespace pargeo;
  using namespace pargeo::dynKdTree;
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state) {
    euclideanMst<dim>(P);
  }
}

static void hull_2d(benchmark::State& state) {
  using namespace pargeo;
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state) {
    pargeo::hull2d::divideConquer::compute(parlay::make_slice(P));
  }
}

static void hull_3d(benchmark::State& state) {
  using namespace pargeo;
  static const int dim = 3;
  auto P = data3d();

  for (auto _ : state) {
    pargeo::hull3d::pseudo::compute(parlay::make_slice(P));
  }
}

static void seb_2d(benchmark::State& state) {
  using namespace pargeo;
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state) {
    pargeo::seb::sampling::compute<dim>(parlay::make_slice(P));
  }
}

static void seb_5d(benchmark::State& state) {
  using namespace pargeo;
  static const int dim = 5;
  auto P = data5d();

  for (auto _ : state) {
    pargeo::seb::sampling::compute<dim>(parlay::make_slice(P));
  }
}

static void knnGraph_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state)
    parlay::sequence<pargeo::dirEdge> E = pargeo::knnGraph(P, 1);
}

static void delaunayGraph_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::delaunayGraph(P);
}

static void gabrielGraph_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::gabrielGraph(P);
}

static void betaSkeleton_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::betaSkeleton(P, 1);
}

static void spanner_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();

  for (auto _ : state)
    parlay::sequence<pargeo::edge> E = pargeo::spanner(P, 40);
}

static void closestPair_2d(benchmark::State& state) {
  static const int dim = 2;
  auto P = data2d();
  for (auto _ : state) pargeo::closestPair(P);
}

static void closestPair_5d(benchmark::State& state) {
  static const int dim = 5;
  auto P = data5d();
  for (auto _ : state) pargeo::closestPair(P);
}

static void mortonSort_2d(benchmark::State& state) {
  auto P = data2d();
  for (auto _ : state) pargeo::zorderSort2d(P);
}

static void mortonSort_3d(benchmark::State& state) {
  auto P = data3d();
  for (auto _ : state) pargeo::zorderSort3d(P);
}


BENCHMARK(kdTree_build_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_knn_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_rangeSearch_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_orthRangeSearch_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_bccp_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_wspd_s2_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_build_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_knn_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_rangeSearch_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_orthRangeSearch_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_bccp_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(kdTree_wspd_s2_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(dynamic_kdTree_build_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(dynamic_kdTree_insert_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(dynamic_kdTree_erase_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(emst_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(hull_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(hull_3d)->Unit(benchmark::kMillisecond);
BENCHMARK(seb_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(seb_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(closestPair_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(closestPair_5d)->Unit(benchmark::kMillisecond);
BENCHMARK(mortonSort_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(mortonSort_3d)->Unit(benchmark::kMillisecond);
BENCHMARK(knnGraph_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(delaunayGraph_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(gabrielGraph_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(betaSkeleton_2d)->Unit(benchmark::kMillisecond);
BENCHMARK(spanner_2d)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
