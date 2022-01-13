#include "pargeo/pointIO.h"
#include "kdTree/kdTree.h"
#include "gtest/gtest.h"

template <int dim>
inline void testKdTree(parlay::sequence<pargeo::point<dim>>& P) {
  using namespace pargeo;
  using nodeT = kdTree::node<dim, point<dim>>;

  nodeT* tree1 = kdTree::build<dim, point<dim>>(P, true);
  nodeT* tree2 = kdTree::build<dim, point<dim>>(P, false);
  nodeT* tree3 = kdTree::build<dim, point<dim>>(P, true, 1);
  nodeT* tree4 = kdTree::build<dim, point<dim>>(P, false, 1);

  std::function<size_t(nodeT*)> checkSum =
    [&](nodeT* node)->size_t {
      if (!node->isLeaf()) {
	size_t lSize = checkSum(node->L());
	size_t rSize = checkSum(node->R());

	// Check if node sizes are consistent
	EXPECT_EQ(lSize + rSize, node->size());

	// Check if each point is in the bounding box
	for (size_t i = 0; i < node->size(); ++ i) {
	  auto p = *(node->at(i));
	  for (size_t d = 0; d < node->dim; ++ d) {
	    EXPECT_TRUE(p[d] <= node->getMax(d));
	    EXPECT_TRUE(p[d] >= node->getMin(d));
	  }
	}

	// Check if box sizes are consistent
	for (size_t d = 0; d < node->dim; ++ d) {
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

TEST(kdTree_structure, testSerial2d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/2d_100.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<2>> P =
    readPointsFromFile<pargeo::point<2>>(filePath);

  testKdTree<2>(P);
}

TEST(kdTree_structure, testSerial3d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/3d_40.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<3>> P =
    readPointsFromFile<pargeo::point<3>>(filePath);

  testKdTree<3>(P);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
