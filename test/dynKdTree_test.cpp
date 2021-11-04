#include <memory>

#include "dataset/uniform.h"
#include "pargeo/pointIO.h"
#include "pargeo/dynKdTree.h"
#include "gtest/gtest.h"

template <int dim>
inline void testKdTree(parlay::sequence<pargeo::point<dim>>& P) {
  using namespace pargeo::dynKdTree;
  using namespace pargeo;

  using nodeT = node<dim, point<dim>>;

  nodeT *tree1 = new nodeT(P); // todo rename the node

  delete tree1;
}

template <int dim>
parlay::sequence<pargeo::point<dim>> data() {
  return pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(200, 0, 1.0);
}

TEST(kdTree_structure, testSerial2d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  static const int dim = 2;

  parlay::sequence<pargeo::point<dim>> P = data<dim>();

  testKdTree<dim>(P);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
