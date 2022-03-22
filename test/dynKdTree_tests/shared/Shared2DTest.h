#ifndef TEST_SHARED2DTEST_H
#define TEST_SHARED2DTEST_H

#include "BasicStructure.h"
#include <gtest/gtest.h>
//#include "common/geometryIO.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"

#include <algorithm>
//#include <kdtree/shared/box.h>
#include "dynamicKdTree/shared/box.h"

typedef pargeo::point<2> pointT;
template <typename Tree>
class Shared2DTest : public BasicStructure2D<Tree> {
 public:
  static const int DIM = 2;
};

TYPED_TEST_SUITE_P(Shared2DTest);

TYPED_TEST_P(Shared2DTest, Verify) {
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  ASSERT_TRUE(tree.verify());
  auto points = this->RESOURCES_1000();
  for (auto& p : points) {
    ASSERT_TRUE(tree.contains(p));
  }

  double point_buf[2];
  for (int i = 0; i < 100; i++) {
    point_buf[0] = (double)i;
    point_buf[1] = (double)i;
    pargeo::point<this->DIM> p(point_buf);
    ASSERT_FALSE(tree.contains(p));
  }
}

TYPED_TEST_P(Shared2DTest, SimpleDelete) {
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();

  // make sure it's built correctly
  for (auto& p : points)
    ASSERT_TRUE(tree.contains(p));

  // Individual Deletes
  for (size_t i = 0; i < points.size() / 2; i += 10) {
    ASSERT_TRUE(tree.erase(points[i]));
  }

  // Bulk Delete
  const auto bulk_delete_size = 25;
  auto start_it = points.begin() + points.size() / 2;
  auto end_it = start_it + bulk_delete_size;
  tree.template erase<false>(parlay::slice(start_it, end_it));

  // Verify
  for (size_t i = 0; i < points.size() / 2; i++) {
    if (i % 10 == 0)
      ASSERT_FALSE(tree.contains(points[i]));
    else
      ASSERT_TRUE(tree.contains(points[i]));
  }
  for (auto i = start_it; i < end_it; i++) {
    ASSERT_FALSE(tree.contains(*i));
  }
  for (auto i = end_it; i < points.end(); i++) {
    ASSERT_TRUE(tree.contains(*i));
  }
}

TYPED_TEST_P(Shared2DTest, SerialDelete) {
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();
  auto to_remove = KEEP_EVEN(points);

  // bulk delete the points
  tree.template erase<false>(to_remove);

  // verify the deletion
  for (size_t i = 0; i < points.size(); i++) {
    if (i % 2 == 0)
      ASSERT_FALSE(tree.contains(points[i]));
    else
      ASSERT_TRUE(tree.contains(points[i]));
  }
}

TYPED_TEST_P(Shared2DTest, BulkDelete) {
  auto tree = this->CONSTRUCT_RESOURCES_1000();
  auto points = this->RESOURCES_1000();
  auto to_remove = KEEP_EVEN(points);

  // bulk delete the points
  tree.bulk_erase(to_remove);

  // verify the deletion
  for (size_t i = 0; i < points.size(); i++) {
    if (i % 2 == 0)
      ASSERT_FALSE(tree.contains(points[i]));
    else
      ASSERT_TRUE(tree.contains(points[i]));
  }
}

TYPED_TEST_P(Shared2DTest, BulkInsert) {
  auto tree = this->CONSTRUCT_RESOURCES_1000_EVEN();
  auto points = this->RESOURCES_1000();
  auto to_insert = KEEP_ODD(points);

  // make sure even points aren't in tree
  for (auto& pt : to_insert) {
    ASSERT_FALSE(tree.contains(pt));
  }

  // bulk delete the points
  tree.insert(to_insert.cut(0, to_insert.size()));

  // verify the insertion
  int i = 0;
  for (auto& pt : points) {
    ASSERT_TRUE(tree.contains(pt)) << "POINT #" << i;
    i++;
  }
}

REGISTER_TYPED_TEST_SUITE_P(
    Shared2DTest, Verify, SimpleDelete, SerialDelete, BulkDelete, BulkInsert);

#endif  // TEST_SHARED2DTEST_H