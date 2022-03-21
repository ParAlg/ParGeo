#ifndef TEST_BINARYHEAPLAYOUT_BHLSTRUCTURE2D_H
#define TEST_BINARYHEAPLAYOUT_BHLSTRUCTURE2D_H

#include <gtest/gtest.h>
//#include "common/geometryIO.h"
#include "pargeo/pointIO.h"

// #include <kdtree/binary-heap-layout/bhlkdtree.h>
// #include <kdtree/shared/macro.h>
#include "dynamicKdTree/binary-heap-layout/bhlkdtree.h"
#include "dynamicKdTree/shared/macro.h"
#include "../shared/BasicStructure.h"

// #include "../shared/BasicStructure.h"

template <typename Tree>
class BHL2DStructureTest : public BasicStructure2D<Tree> {
 public:
  static const int DIM = 2;
};

TYPED_TEST_SUITE_P(BHL2DStructureTest);

TYPED_TEST_P(BHL2DStructureTest, LayoutSize2) {
  // create the tree
  auto tree = this->CONSTRUCT_2D_SIZE_2();
  auto root = tree.root();  // get the root to verify the layout

  // Check node types
  ASSERT_FALSE(root[0].isLeaf());
  ASSERT_TRUE(root[1].isLeaf());
  ASSERT_TRUE(root[2].isLeaf());
  // Check point counts
  ASSERT_EQ(root[0].countPoints(), 2);
  ASSERT_EQ(root[1].countPoints(), 1);
  ASSERT_EQ(root[2].countPoints(), 1);
  // Check memory values
  ASSERT_EQ(root[0].getSplitDimension(), 0);
  //#ifdef USE_MEDIAN_SELECTION
  ASSERT_EQ(root[0].getSplitValue(), 1.0);
  //#else
  // ASSERT_EQ(root[0].getSplitValue(), 0.5);
  //#endif
  ASSERT_EQ(root[1].getValues().size(), 1);
  ASSERT_EQ(root[1].getValues()[0], this->POINT_ARR_2[0]);
  ASSERT_EQ(tree.getNodeValueIdx(&root[1]), std::make_pair(0, 1));
  ASSERT_EQ(root[2].getValues().size(), 1);
  ASSERT_EQ(root[2].getValues()[0], this->POINT_ARR_2[1]);
  ASSERT_EQ(tree.getNodeValueIdx(&root[2]), std::make_pair(1, 2));
  // Check memory layout
  ASSERT_EQ(root->getLeft(), root + 1);
  ASSERT_EQ(root->getRight(), root + 2);
}

TYPED_TEST_P(BHL2DStructureTest, LayoutSize8) {
  // create the tree
  auto tree = this->CONSTRUCT_2D_SIZE_8();
  auto root = tree.root();  // get the root to verify the layout

  // Check node types

  /*               0
   *           1       2
   *         3   4   5   6
   *        7 8 9 A B C D E
   */
  ASSERT_FALSE(root[0].isLeaf());

  ASSERT_FALSE(root[1].isLeaf());
  ASSERT_FALSE(root[2].isLeaf());

  for (int i = 3; i < 7; i++) {
    ASSERT_FALSE(root[i].isLeaf());
  }

  for (int i = 7; i < 15; i++) {
    ASSERT_TRUE(root[i].isLeaf());
  }

  // Check point counts
  ASSERT_EQ(root[0].countPoints(), 8);

  ASSERT_EQ(root[1].countPoints(), 4);
  ASSERT_EQ(root[2].countPoints(), 4);

  for (int i = 3; i < 7; i++) {
    ASSERT_EQ(root[i].countPoints(), 2);
  }

  for (int i = 7; i < 15; i++) {
    ASSERT_EQ(root[i].countPoints(), 1);
  }

  // Check memory values - serial case always uses median selection
  //#ifdef USE_MEDIAN_SELECTION
  double split_values[7] = {4, 2, 6, 1, 3, 5, 7};
  //#else
  // double split_values[7] = {3.5, 1.5, 5.5, 0.5, 2.5, 4.5, 6.5};
  //#endif

  ASSERT_EQ(root[0].getSplitDimension(), 0);
  ASSERT_EQ(root[0].getSplitValue(), split_values[0]);

  ASSERT_EQ(root[1].getSplitDimension(), 1);
  ASSERT_EQ(root[1].getSplitValue(), split_values[1]);
  ASSERT_EQ(root[2].getSplitDimension(), 1);
  ASSERT_EQ(root[2].getSplitValue(), split_values[2]);

  ASSERT_EQ(root[3].getSplitDimension(), 0);
  ASSERT_EQ(root[3].getSplitValue(), split_values[3]);
  ASSERT_EQ(root[4].getSplitDimension(), 0);
  ASSERT_EQ(root[4].getSplitValue(), split_values[4]);
  ASSERT_EQ(root[5].getSplitDimension(), 0);
  ASSERT_EQ(root[5].getSplitValue(), split_values[5]);
  ASSERT_EQ(root[6].getSplitDimension(), 0);
  ASSERT_EQ(root[6].getSplitValue(), split_values[6]);

  for (int i = 0; i < 8; i++) {
    ASSERT_EQ(root[i + 7].getValues().size(), 1);
    ASSERT_EQ(root[i + 7].getValues()[0], this->POINT_ARR_8[i]);
    ASSERT_EQ(tree.getNodeValueIdx(&root[i + 7]), std::make_pair(i, i + 1));
  }

  // Check memory layout
  for (int i = 0; i < 7; i++) {
    ASSERT_EQ(root[i].getLeft(), root + 2 * i + 1);
    ASSERT_EQ(root[i].getRight(), root + 2 * i + 2);
  }
}

REGISTER_TYPED_TEST_SUITE_P(BHL2DStructureTest, LayoutSize2, LayoutSize8);

#endif  // TEST_BINARYHEAPLAYOUT_BHLSTRUCTURE2D_H
