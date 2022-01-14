#include "dataset/uniform.h"
#include "pargeo/pointIO.h"
#include "kdTree/kdTree.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::point<2>> smallData()
{
  parlay::sequence<pargeo::point<2>> P(7);

  P[0][0] = 0.5;
  P[0][1] = 3.5;
  P[1][0] = 0.5;
  P[1][1] = 2.5;
  P[2][0] = 1;
  P[2][1] = 3;
  P[3][0] = 1.5;
  P[3][1] = 2.5;
  P[4][0] = 3.5;
  P[4][1] = 0.5;
  P[5][0] = 3.5;
  P[5][1] = 0.5;
  P[6][0] = 2.5;
  P[6][1] = 0.5;

  return P;
}

template <int dim>
inline void testKdTree(parlay::sequence<pargeo::point<dim>> &P)
{
  using namespace pargeo;
  using nodeT = kdTree::node<dim, point<dim>>;

  nodeT *tree1 = kdTree::build<dim, point<dim>>(P, true);
  nodeT *tree2 = kdTree::build<dim, point<dim>>(P, false);
  nodeT *tree3 = kdTree::build<dim, point<dim>>(P, true, 1);
  nodeT *tree4 = kdTree::build<dim, point<dim>>(P, false, 1);

  std::function<size_t(nodeT *)> checkSum =
      [&](nodeT *node) -> size_t
  {
    if (!node->isLeaf())
    {
      size_t lSize = checkSum(node->L());
      size_t rSize = checkSum(node->R());

      // Check if node sizes are consistent
      EXPECT_EQ(lSize + rSize, node->size());

      // Check if each point is in the bounding box
      for (size_t i = 0; i < node->size(); ++i)
      {
        auto p = *(node->at(i));
        for (size_t d = 0; d < node->dim; ++d)
        {
          EXPECT_TRUE(p[d] <= node->getMax(d));
          EXPECT_TRUE(p[d] >= node->getMin(d));
        }
      }

      // Check if box sizes are consistent
      for (size_t d = 0; d < node->dim; ++d)
      {
        EXPECT_FLOAT_EQ(std::max(node->L()->getMax(d),
                                 node->R()->getMax(d)),
                        node->getMax(d));
      }
    }
    return node->size();
  };

  checkSum(tree1);
  checkSum(tree2);
  checkSum(tree3);
  checkSum(tree4);

  kdTree::del(tree1);
  kdTree::del(tree2);
  kdTree::del(tree3);
  kdTree::del(tree4);
}

TEST(kdTree_structure, testSerial2d)
{
  static const int dim = 2;
  auto P = pargeo::uniformInPolyPoints<dim>(100, 0, 1.0);
  testKdTree<dim>(P);
}

TEST(kdTree_structure, testSerial3d)
{
  static const int dim = 5;
  auto P = pargeo::uniformInPolyPoints<dim>(100, 0, 1.0);
  testKdTree<dim>(P);
}

TEST(kdTree_method, knnSearch)
{
  using namespace pargeo;
  using namespace parlay;
  auto P = smallData();

  auto tree = kdTree::build<2, point<2>>(P, false, 1);

  sequence<size_t> nn = kdTree::batchKnn(P, 2, tree, true);

  EXPECT_EQ(P.size(), 7);
  EXPECT_EQ(nn.size(), 14);
  EXPECT_EQ(nn[0], 0);
  EXPECT_EQ(nn[2], 1);
  EXPECT_EQ(nn[4], 2);
  EXPECT_EQ(nn[6], 3);
  // 4th and 5th points are the same -- NN ordering may be arbitrary
  // -> test omitted
  EXPECT_EQ(nn[12], 6);

  EXPECT_FLOAT_EQ(P[0].dist(P[nn[1]]), sqrt(2) / 2);
  EXPECT_FLOAT_EQ(P[1].dist(P[nn[3]]), sqrt(2) / 2);
  EXPECT_FLOAT_EQ(P[2].dist(P[nn[5]]), sqrt(2) / 2);
  EXPECT_FLOAT_EQ(P[3].dist(P[nn[7]]), sqrt(2) / 2);
  EXPECT_FLOAT_EQ(P[4].dist(P[nn[9]]), 0);
  EXPECT_FLOAT_EQ(P[5].dist(P[nn[11]]), 0);
  EXPECT_FLOAT_EQ(P[6].dist(P[nn[13]]), 1);

  kdTree::del(tree);
}

TEST(kdTree_method, rangeSearch)
{
  using namespace pargeo;
  using namespace parlay;
  auto P = smallData();

  auto tree = kdTree::build<2, point<2>>(P, false, 1);

  {
    sequence<point<2> *> nbrs = kdTree::rangeSearch(tree, P[0], sqrt(2) / 2);
    EXPECT_EQ(nbrs.size(), 2);
  }

  {
    sequence<point<2> *> nbrs = kdTree::rangeSearch(tree, P[0], 1);
    EXPECT_EQ(nbrs.size(), 3);
  }

  {
    sequence<point<2> *> nbrs = kdTree::rangeSearch(tree, P[0], sqrt(2));
    EXPECT_EQ(nbrs.size(), 4);
  }

  {
    sequence<point<2> *> nbrs = kdTree::rangeSearch(tree, P[0], sqrt(2) * 3);
    EXPECT_EQ(nbrs.size(), 7);
  }

  kdTree::del(tree);
}

TEST(kdTree_method, orthogonalRangeSearch)
{
  using namespace pargeo;
  using namespace parlay;
  auto P = smallData();

  auto tree = kdTree::build<2, point<2>>(P, false, 1);

  {
    sequence<point<2> *> nbrs = kdTree::orthogonalRangeSearch(tree, P[0], 0.5);
    EXPECT_EQ(nbrs.size(), 2);
  }

  {
    sequence<point<2> *> nbrs = kdTree::orthogonalRangeSearch(tree, P[0], 1);
    EXPECT_EQ(nbrs.size(), 4);
  }

  {
    sequence<point<2> *> nbrs = kdTree::orthogonalRangeSearch(tree, P[0], sqrt(2));
    EXPECT_EQ(nbrs.size(), 4);
  }

  {
    sequence<point<2> *> nbrs = kdTree::orthogonalRangeSearch(tree, P[0], 3);
    EXPECT_EQ(nbrs.size(), 7);
  }

  kdTree::del(tree);
}

TEST(kdTree_method, bichromaticClosestPair)
{
  using namespace pargeo;
  using namespace parlay;

  auto P = smallData();

  auto tree1 = kdTree::build<2, point<2>>(P.cut(0, 4), false, 1);

  auto tree2 = kdTree::build<2, point<2>>(P.cut(4, 7), false, 1);

  auto result = kdTree::bichromaticClosestPair(tree1, tree2);

  EXPECT_FLOAT_EQ(std::get<0>(result)->dist(*std::get<1>(result)), sqrt(5));

  EXPECT_FLOAT_EQ(std::get<2>(result), sqrt(5));
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
