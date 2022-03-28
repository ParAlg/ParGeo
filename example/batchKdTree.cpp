#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "batchKdtree/shared/geometry.h"
#include "pargeo/parseCommandLine.h"
#include "pargeo/pointIO.h"
#include "pargeo/getTime.h"

#include "batchKdtree/cache-oblivious/cokdtree.h"
#include "batchKdtree/binary-heap-layout/bhlkdtree.h"
#include "batchKdtree/log-tree/logtree.h"
#include "batchKdtree/shared/dual.h"

// using namespace benchIO;
using namespace pargeo::batchKdTree;

using coord = double;

enum TestType {
  INVALID = 0,
  CONSTRUCTION,
  RANGE_QUERY,
  CONTAINS,
  KNN,
  KNN2,
  KNN3,
  DUAL_KNN,
  DELETE,
  INSERT,
  INSERT_DELETE
};
bool is_knn(const TestType& t) {
  return (t == KNN) || (t == KNN2) || (t == KNN3) || (t == DUAL_KNN);
}

static std::ostream& operator<<(std::ostream& os, const TestType& t) {
  switch (t) {
    case INVALID:
      os << "INVALID";
      break;
    case CONSTRUCTION:
      os << "CONSTRUCTION";
      break;
    case RANGE_QUERY:
      os << "RANGE_QUERY";
      break;
    case CONTAINS:
      os << "CONTAINS";
      break;
    case KNN:
      os << "KNN";
      break;
    case KNN2:
      os << "KNN2";
      break;
    case KNN3:
      os << "KNN3";
      break;
    case DUAL_KNN:
      os << "DUAL_KNN";
      break;
    case DELETE:
      os << "DELETE";
      break;
    case INSERT:
      os << "INSERT";
      break;
    case INSERT_DELETE:
      os << "INSERT + DELETE";
      break;
  }
  return os;
}

struct TestOptions {
  bool run_serial;
  bool run_parallel;
  bool run_co;
  bool run_bhl;
  bool run_log;
  int rounds;
  int k;           // for knn
  int k_type;      // for knn
  int percentage;  // for delete/construction
  TestType type;
  const char* out_file;

  TestOptions()
      : run_serial(false),
        run_parallel(false),
        run_co(false),
        run_bhl(false),
        run_log(false),
        rounds(0),
        k(3),
        k_type(0),
        percentage(-1),
        type(INVALID),
        out_file(nullptr) {}
};

template <int dim, class Tree>
void timeConstruction(parlay::sequence<point<dim>> const& P, int rounds, int percentage) {
  // P is our local copy, we can move into tree
  if (percentage == -1) percentage = 100;  // default
  auto P_slice = P.cut(0, (P.size() * percentage) / 100);
  // auto construct_size = (P.size() * percentage) / 100;
  // parlay::sequence<point<dim>> P_construct;
  // P_construct.assign(P.begin(), P.begin() + construct_size);
  pargeo::timer t;
  for (int i = 0; i < rounds; ++i) {
    t.start();
    Tree tree(P_slice);
    auto build_time = t.get_next();
    // std::cout << "  -> build-time = " << build_time << std::endl;
    if (i > 0) std::cout << ",  ";
    std::cout << build_time;
  }
  t.stop();
}

template <int dim, class Tree>
void timeContainsQuery(parlay::sequence<point<dim>> const& P, int rounds) {
  constexpr int rep = 10000;
  double x[dim];
  point<dim> pnot[rep];

  for (int j = 0; j < rep; j++) {
    for (int i = 0; i < dim; i++) {
      x[i] = i;
    }
    pnot[j] = point<dim>(x);
  }

  pargeo::timer t;
  for (int i = 0; i < rounds; ++i) {
    if (i > 0) std::cout << ",  ";
    Tree tree(P);

    t.start();
    for (int j = 0; j < rep; j++) {
      tree.contains(pnot[j]);
    }
    auto query_time_1 = t.get_next();
    // std::cout << "  -> " << rep << "x query-time (-) = " << query_time_1 << std::endl;
    std::cout << "[" << query_time_1 << ", ";
    t.reset();

    t.start();
    for (int j = 0; j < rep; j++) {
      tree.contains(P[j]);
    }
    auto query_time_2 = t.get_next();
    // std::cout << "  -> " << rep << "x query-time (+) = " << query_time_2 << std::endl;
    std::cout << query_time_2 << "]";
  }
  t.stop();
}

template <int dim, class Tree>
void timeRangeQuery(parlay::sequence<point<dim>> const& P, int rounds) {
  constexpr int rep = 10;

  // compute full bounding box
  point<dim> pMin, pMax;
  pMin = P[0];
  pMax = P[0];
  for (const auto& pt : P) {
    pMin.minCoords(pt);
    pMax.maxCoords(pt);
  }
  /*auto printPoint = [](const point<dim>& p) {
    std::stringstream ss;
    ss << "(";
    for (int i = 0; i < dim; i++) {
      if (i > 0) ss << ", ";
      ss << p.coordinate(i);
    }
    ss << ")";
    return ss.str();
  };
  std::cout << "Min: " << printPoint(pMin) << std::endl;
  std::cout << "Max: " << printPoint(pMax) << std::endl;*/

  // scale in by 10x on each dimension -> 1% for 2d
  double x[dim];
  for (int i = 0; i < dim; i++)
    x[i] = pMin[i] / 10;
  auto qMin = point<dim>(x);
  for (int i = 0; i < dim; i++)
    x[i] = pMax[i] / 10;
  auto qMax = point<dim>(x);

  pargeo::timer t;
  for (int i = 0; i < rounds; ++i) {
    Tree tree(P);

    parlay::sequence<point<dim>> p;
    t.start();
    for (int j = 0; j < rep; j++) {
      p = tree.orthogonalQuery(qMin, qMax);
    }
    auto time = t.get_next();
    // std::cout << "  -> " << rep << "x query-time (" << p.size() << " points) = " <<time
    //<< std::endl;
    if (i > 0) std::cout << ", ";
    std::cout << time;
  }
  t.stop();
}

template <int dim, class Tree>
void timeKnn(parlay::sequence<point<dim>>& P, int rounds, int k, int k_type, const TestType& type) {
  // constexpr int rep = 1;

  for (int i = 0; i < rounds; ++i) {
    if (i > 0) std::cout << ", ";
    Tree tree(P);

    pargeo::timer t;
    t.start();
    // for (int j = 0; j < rep; j++) {
    if (type == KNN) {
      if (k_type == 0) tree.template knn<false, false>(P, k);
      if (k_type == 1) tree.template knn<false, true>(P, k);
      if (k_type == 2) tree.template knn<true, false>(P, k);
      if (k_type == 3) tree.template knn<true, true>(P, k);
    } else if (type == KNN2) {
      if (k_type == 0) tree.template knn2<false, false>(P, k);
      if (k_type == 1) tree.template knn2<false, true>(P, k);
      if (k_type == 2) tree.template knn2<true, false>(P, k);
      if (k_type == 3) tree.template knn2<true, true>(P, k);
    } else if (type == KNN3) {
      if (k_type == 0) tree.template knn3<false, false>(P, k);
      if (k_type == 1) tree.template knn3<false, true>(P, k);
      if (k_type == 2) tree.template knn3<true, false>(P, k);
      if (k_type == 3) tree.template knn3<true, true>(P, k);
    } else {
      dualKnn(P, tree, k);
    }
    //}
    auto time = t.get_next();
    std::cout << time;
    // std::cout << "  -> " << rep << "x query-time (k=" << k << ") = " << time << std::endl;
    t.stop();
  }
}

template <int dim, class Tree>
void timeInsert(parlay::sequence<point<dim>> const& P, int rounds) {
  int log2size = (int)std::ceil(std::log2(P.size()));
  constexpr int divisor = 10;
  int div_size = P.size() / divisor;

  for (int i = 0; i < rounds; ++i) {
    if (i > 0) std::cout << ", ";
    Tree tree(log2size);

    pargeo::timer t;
    t.start();
    for (int j = 0; j < divisor; j++) {
      tree.insert(P.cut(j * div_size, (j + 1) * div_size));
    }
    auto time = t.get_next();
    std::cout << time;
    t.stop();
  }
}

template <int dim, class Tree>
void timeDelete(parlay::sequence<point<dim>> const& P, int rounds, int percentage) {
  if (percentage == -1) percentage = 10;  // default

  constexpr int divisor = 10;
  if (percentage % divisor != 0) throw std::runtime_error("delete setup error!");
  int div_size = (P.size() * percentage) / 100;
  parlay::sequence<point<dim>> to_delete;
  to_delete.reserve(div_size);
  for (size_t i = 0; i < P.size(); i++) {
    if (i % divisor < (size_t)(percentage / divisor)) {
      to_delete.push_back(P[i]);
    }
  }

  for (int i = 0; i < rounds; ++i) {
    if (i > 0) std::cout << ", ";
    Tree tree(P);

    pargeo::timer t;
    t.start();
    tree.bulk_erase(to_delete);
    auto time = t.get_next();
#ifdef ERASE_SEARCH_TIMES
    std::cout << "TOTAL SEARCH TIME: " << tree.total_search_time << " s" << std::endl;
    std::cout << "TOTAL remove TIME: " << tree.total_bbox_time << " s" << std::endl;
    std::cout << "TOTAL LEAF TIME: " << tree.total_leaf_time << " s" << std::endl;
#endif
    std::cout << time;
    t.stop();
  }
}

// /*template <int dim, class Tree>
// void timeInsertDelete(parlay::sequence<point<dim>> const& P, int rounds) {
//   int log2size = (int)std::ceil(std::log2(P.size()));

//   for (int i = 0; i < rounds; ++i) {
//     if (i > 0) std::cout << ", ";
//     Tree tree(log2size);

//     timer t;
//     t.start();

//     t.stop();
//   }
// }*/

template <int dim, class Tree>
void dispatchTest(parlay::sequence<point<dim>>& P, const TestOptions& test_options) {
  std::cout << "[";
  switch (test_options.type) {
    case CONSTRUCTION: {
      timeConstruction<dim, Tree>(P, test_options.rounds, test_options.percentage);
      break;
    }
    case CONTAINS: {
      timeContainsQuery<dim, Tree>(P, test_options.rounds);
      break;
    }
    case RANGE_QUERY: {
      timeRangeQuery<dim, Tree>(P, test_options.rounds);
      break;
    }
    case DUAL_KNN:
    case KNN3:
    case KNN2:
    case KNN: {
      timeKnn<dim, Tree>(
          P, test_options.rounds, test_options.k, test_options.k_type, test_options.type);
      break;
    }
    case DELETE: {
      timeDelete<dim, Tree>(P, test_options.rounds, test_options.percentage);
      break;
    }
    case INSERT: {
      timeInsert<dim, Tree>(P, test_options.rounds);
      break;
    }
    // case INSERT_DELETE: {
    // timeInsertDelete<dim, Tree>(P, test_options.rounds);
    // break;
    //}
    default: {
      throw std::runtime_error("Unsupported test type!");
      break;
    }
  }
  std::cout << "]" << std::endl;
}

static bool UseCO(const TestType& t) {
  return (t != INSERT) && (t != INSERT_DELETE) && (t != KNN2) && (t != KNN3);
}
static bool UseBHL(const TestType& t) { return (t != KNN2) && (t != KNN3); }

template <int dim, bool parallel, bool coarsen>
void timeTrees(parlay::sequence<point<dim>>& P, const TestOptions& test_options) {
  typedef CO_KdTree<dim, point<dim>, parallel, coarsen> COTree_t;
  typedef BHL_KdTree<dim, point<dim>, parallel, coarsen> BHLTree_t;
  constexpr int NUM_TREES = 21;
  constexpr int BUFFER_LOG2_SIZE = 7;
  typedef LogTree<NUM_TREES, BUFFER_LOG2_SIZE, dim, point<dim>, parallel, coarsen> LogTree_t;

  std::string thread_str =
      (parallel ? (" (" + std::to_string(parlay::num_workers()) + " threads)") : (" (serial)"));
  if (test_options.run_co && UseCO(test_options.type)) {
    std::cout << " - Timing Cache-Oblivious Tree" << thread_str << std::endl;
    dispatchTest<dim, COTree_t>(P, test_options);
  }

  if (test_options.run_bhl && UseBHL(test_options.type)) {
    std::cout << " - Timing Binary Heap Layout Tree" << thread_str << std::endl;
    dispatchTest<dim, BHLTree_t>(P, test_options);
  }

  if (test_options.run_log) {
    std::cout << " - Timing Log Tree" << thread_str << std::endl;

    // print tree mask for debugging
#ifndef NDEBUG
    LogTree_t logtree(P);
    std::cout << " -> tree_mask = " << logtree.getTreeMask() << std::endl;
#endif

    // time it
    dispatchTest<dim, LogTree_t>(P, test_options);
  }
}

template <int dim>
inline void runTests(const char* iFile, const TestOptions& test_options) {
  // parlay::sequence<point<dim>> Points = readPointsFromFile<point<dim>>(iFile);
  parlay::sequence<point<dim>> Points = pargeo::pointIO::readPointsFromFile<point<dim>>(iFile);
  std::cout << "Timing " << test_options.type << " (dim = " << dim
            << ", #points = " << Points.size()
            << ((is_knn(test_options.type)) ? (", k=" + std::to_string(test_options.k)) : "") << ")"
            << std::endl;

  constexpr bool coarsen = false;
  if (test_options.run_serial) {
    // std::cout << "Serial Timings ----------------" << std::endl;
    timeTrees<dim, false, coarsen>(Points, test_options);
  }
  if (test_options.run_parallel) {
    // std::cout << "Parallel Timings (#threads = " << parlay::num_workers() << ") ----------------"
    //<< std::endl;
    timeTrees<dim, true, coarsen>(Points, test_options);
  }
}

TestOptions parseCmdLine(pargeo::commandLine& P) {
  // parse test options
  TestOptions test_options;
  test_options.out_file = P.getOptionValue("-o");
  test_options.rounds = P.getOptionIntValue("-r", 1);
  test_options.k = P.getOptionIntValue("-k", 3);
  test_options.k_type = P.getOptionIntValue("--knn_type", 0);
  test_options.percentage = P.getOptionIntValue("--percentage", -1);

  // run type
  test_options.run_serial = P.getOption("--serial");
  test_options.run_parallel = P.getOption("--parallel");
  if (!(test_options.run_serial || test_options.run_parallel)) {
    test_options.run_serial = test_options.run_parallel = true;  // neither set, so default to both
  }

  test_options.run_co = P.getOption("--co");
  test_options.run_bhl = P.getOption("--bhl");
  test_options.run_log = P.getOption("--log");
  if (!(test_options.run_co || test_options.run_bhl || test_options.run_log)) {
    // none set, so default to all
    test_options.run_co = test_options.run_bhl = test_options.run_log = true;
  }

  auto type_str = P.getOptionValue("--type", "construction");  // default to construction times
  if (type_str == "construction") {
    test_options.type = CONSTRUCTION;
  } else if (type_str == "range") {
    test_options.type = RANGE_QUERY;
  } else if (type_str == "contains") {
    test_options.type = CONTAINS;
  } else if (type_str == "knn") {
    test_options.type = KNN;
  } else if (type_str == "knn2") {
    test_options.type = KNN2;
  } else if (type_str == "knn3") {
    test_options.type = KNN3;
  } else if (type_str == "dual_knn") {
    test_options.type = DUAL_KNN;
  } else if (type_str == "delete") {
    test_options.type = DELETE;
  } else if (type_str == "insert") {
    test_options.type = INSERT;
  } else if (type_str == "insert_delete") {
    test_options.type = INSERT_DELETE;
  } else {
    throw std::runtime_error("Invalid type: " + type_str);
  }

  return test_options;
}

int main(int argc, char* argv[]) {
  pargeo::commandLine P(
      argc,
      argv,
      "[-o <outFile>] [-r <rounds>] [--serial, --parallel] [ --co, --bhl, "
      "--log] [--type "
      "construction|range|contains|knn|knn2|knn3|dual_knn|delete|insert|insert_delete] [-k "
      "<k for knn>] [--knn_type <[0,3]>] [--percentage <percentage for "
      "construction/deletion>] <inFile>");
  char* iFile = P.getArgument(0);
  auto test_options = parseCmdLine(P);

  int dim = pargeo::pointIO::readHeader(iFile);
  if (dim == 2) {
    runTests<2>(iFile, test_options);
  } else if (dim == 3) {
    runTests<3>(iFile, test_options);
  }
}
