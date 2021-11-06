#include <memory>

#include "dataset/uniform.h"
#include "parlay/sequence.h"
#include "pargeo/pointIO.h"
#include "pargeo/dynKdTree.h"
#include "gtest/gtest.h"


template <int dim>
parlay::sequence<pargeo::point<dim>> data() {
  return pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(200, 0, 1.0);
}


template <int dim>
inline void testKdTree(parlay::sequence<pargeo::point<dim>>& P) {
  using namespace pargeo::dynKdTree;

  using nodeT = rootNode<dim, pargeo::point<dim>>;

  std::unique_ptr<nodeT> tree1 = std::unique_ptr<nodeT>(new nodeT(P));

  tree1->insert(P);

  tree1->erase(P);

  EXPECT_TRUE(tree1->check());

  tree1->kNN(P[0], 10);
}


TEST(dynKdTree, treeStructure) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  static const int dim = 2;

  parlay::sequence<pargeo::point<dim>> P = data<dim>();

  testKdTree<dim>(P);
}


TEST(dynKdTree, boundingBox) {
  using namespace pargeo::dynKdTree;
  static const int dim = 2;

  double data[12] = {
    0.1, 0.1,
    0.3, 0.3,
    0.0, 0.0,
    0.4, 0.4,
    0.31, 0.1,
    0.4, 0.2
  };

  parlay::sequence<pargeo::point<dim>> b1;
  b1.emplace_back(data);
  b1.emplace_back(data + 2);

  parlay::sequence<pargeo::point<dim>> b2;
  b2.emplace_back(data + 4);
  b2.emplace_back(data + 6);

  parlay::sequence<pargeo::point<dim>> b3;
  b3.emplace_back(data + 8);
  b3.emplace_back(data + 10);

  auto bb1 = boundingBox<dim>(b1);
  auto bb2 = boundingBox<dim>(b2);
  auto bb3 = boundingBox<dim>(b3);

  ASSERT_FALSE(bb1.contains(b2[0]));
  ASSERT_TRUE(bb2.contains(b1[0]));
  ASSERT_TRUE(bb1.compare(bb2) == boundingBox<dim>::overlap);
  ASSERT_TRUE(bb2.compare(bb1) == boundingBox<dim>::include);
  ASSERT_TRUE(bb1.compare(bb3) == boundingBox<dim>::exclude);
  ASSERT_TRUE(bb3.compare(bb1) == boundingBox<dim>::exclude);
  ASSERT_TRUE(bb2.compare(bb3) == boundingBox<dim>::include);
  ASSERT_TRUE(bb3.compare(bb2) == boundingBox<dim>::overlap);
}


TEST(dynKdTree, kBuffer) {
  using namespace pargeo::dynKdTree;

  auto kBuf = kBuffer<int>(4);

  kBuf.insertK({0.1, 1});
  kBuf.insertK({0.4, 1});
  kBuf.insertK({0.6, 1});
  kBuf.insertK({0.2, 1});
  kBuf.insertK({0.5, 1});
  kBuf.insertK({0.3, 1});

  ASSERT_EQ(kBuf.getK().first, 0.4);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
