#ifndef TEST_QUERYTEST_H
#define TEST_QUERYTEST_H

#include "BasicStructure.h"
#include <gtest/gtest.h>
//#include "common/geometryIO.h"
#include "pargeo/pointIO.h"

#include <algorithm>
// #include <kdtree/shared/box.h>
// #include <kdtree/shared/dual.h>
#include "dynamicKdTree/shared/box.h"
#include "dynamicKdTree/shared/dual.h"

typedef pargeo::point<2> pointT;
template <typename Tree>
class QueryTest : public BasicStructure2D<Tree> {
 public:
  static const int DIM = 2;
};

TYPED_TEST_SUITE_P(QueryTest);

TYPED_TEST_P(QueryTest, BasicRangeQuery) {
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();

  pointT qMin, qMax;

  auto compare = [&](pointT& l, pointT& r) {
    return l[0] < r[0];
  };
  parlay::sort_inplace(points, compare);

  // test 1
  for (int i = 0; i < 4; i++) {
    auto div = (1 << i);
    auto size = points.size() / div;
    for (int b = 0; b < div; b++) {
      // get the slice
      auto slice1 = points.cut(b * size, (b + 1) * size);
      boundingBoxSerial(qMin, qMax, slice1);

      // get the query
      auto res1 = tree.orthogonalQuery(qMin, qMax);

      // validate the query
      ASSERT_EQ(res1.size(), slice1.size());
      parlay::sort_inplace(res1, compare);
      for (size_t j = 0; j < slice1.size(); j++)
        ASSERT_EQ(slice1[j], res1[j]);
    }
  }
}

TYPED_TEST_P(QueryTest, BasicKnn) {
  // construct tree
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();
  ASSERT_EQ(tree.size(), points.size());

  // construct brute force solution
  constexpr int k = 4;
  auto check = knnBuf::bruteforceKnn(points, k);
  ASSERT_EQ(check.size(), k * points.size());

  for (int i = 0; i < 4; i++) {
    decltype(check) res;
    if (i == 0) {
      res = tree.template knn<false, false>(points, k);
    } else if (i == 1) {
      res = tree.template knn<false, true>(points, k);
    } else if (i == 2) {
      res = tree.template knn<true, false>(points, k);
    } else {
      res = tree.template knn<true, true>(points, k);
    }

    // verify result
    ASSERT_EQ(res.size(), k * points.size());
    for (size_t i = 0; i < points.size(); i++) {
      // std::cout << "CHECKING knn FOR " << i << ": " << points[i] << std::endl;
      // sort the results for this point
      auto res_start = res.begin() + i * k;
      auto check_start = check.begin() + i * k;
      auto compare = [&](pointT* l, pointT* r) {
        return l->at(0) < r->at(0);
      };
      std::sort(check_start, check_start + k, compare);
      std::sort(res_start, res_start + k, compare);

      for (int j = 0; j < k; j++) {
        auto res_pt = **(res_start + j);
        auto check_pt = **(check_start + j);
        // std::cout << " - " << j << ": " << res_pt << " vs " << check_pt << std::endl;
        EXPECT_EQ(res_pt, check_pt) << "dist res, dist check =  " << points[i].dist(res_pt) << ", "
                                    << points[i].dist(check_pt);
      }
    }
  }
}

TYPED_TEST_P(QueryTest, DualKnn) {
  constexpr int k = 4;
  // construct tree
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();
  ASSERT_EQ(tree.size(), points.size());

  auto res = dualKnn(points, tree, k);
  ASSERT_EQ(res.size(), k * points.size());
  ASSERT_EQ(points.size(), 1000);

  // construct brute force solution
  // HAVE to do this after dualKnn, because it reorders the query points
  auto check = knnBuf::bruteforceKnn(points, k);
  ASSERT_EQ(check.size(), k * points.size());

  // verify result
  for (size_t i = 0; i < points.size(); i++) {
    // std::cout << "CHECKING knn FOR " << i << ": " << points[i] << std::endl;
    // sort the results for this point
    auto res_start = res.begin() + i * k;
    auto check_start = check.begin() + i * k;
    auto compare = [&](pointT* l, pointT* r) {
      return l->at(0) < r->at(0);
    };
    std::sort(check_start, check_start + k, compare);
    std::sort(res_start, res_start + k, compare);

    for (int j = 0; j < k; j++) {
      auto res_pt = **(res_start + j);
      auto check_pt = **(check_start + j);
      // std::cout << " - " << j << ": " << res_pt << " vs " << check_pt << std::endl;
      EXPECT_EQ(res_pt, check_pt) << i << ". (" << j
                                  << ") dist res, dist check =  " << points[i].dist(res_pt) << ", "
                                  << points[i].dist(check_pt);
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(QueryTest, BasicRangeQuery, BasicKnn, DualKnn);

#endif  // TEST_QUERYTEST_H