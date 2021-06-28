#include <tuple>
#include "gtest/gtest.h"
#include "pargeo/octTree.h"
#include "dataset/uniform.h"

template<class objT, class pointT>
class octTreeTest : public pargeo::octTree<objT, pointT> {
  using T = pargeo::octTree<objT, pointT>;
  using floatT = typename T::floatT;
  using node = typename T::node;

  static constexpr floatT eps = 0.0001;

  size_t checkSumHelper(size_t l, size_t nodeIdx) {
    node* cur = &T::nodes[l][nodeIdx];

    if (T::lastLevel(l) || cur->isLeaf()) return cur->size();

    size_t total = 0;
    for (size_t i = 0; i < 8; ++ i) {
      size_t nodeIdxNext = nodeIdx << T::dim;
      nodeIdxNext = nodeIdxNext + i;
      total += checkSumHelper(l + 1, nodeIdxNext);
    }

    EXPECT_EQ(total, cur->size());

    return cur->size();
  }

  void checkElemHelper(size_t l, size_t nodeIdx) {
    node* cur = &T::nodes[l][nodeIdx];

    if (T::lastLevel(l) || cur->isLeaf()) return;

    auto pBox = T::getBox(l, nodeIdx);
    objT pMin = std::get<0>(pBox);
    objT pMax = std::get<1>(pBox);

    for (size_t i = cur->s; i < cur->e; ++ i) {
      objT p = T::items[i].p;
      for (int d = 0; d < 3; ++ d) {
	EXPECT_TRUE(p[d] + eps >= pMin[d]);
	EXPECT_TRUE(p[d] <= pMax[d] + eps);
	if (!(p[d] + eps >= pMin[d]) ||
	    !(p[d] <= pMax[d] + eps)) {
	  std::cout << l << "\n";
	  std::cout << pMin << ", " << pMax << "\n";
	  std::cout << p << "\n";
	  T::printBits(nodeIdx);
	  T::printBits(T::items[i].id);
	  abort();
	}
      }
    }

    for (size_t i = 0; i < 8; ++ i) {
      size_t nodeIdxNext = nodeIdx << T::dim;
      nodeIdxNext = nodeIdxNext + i;
      node* next = &T::nodes[l + 1][nodeIdxNext];
      checkElemHelper(l + 1, nodeIdxNext);
    }
  }

  std::tuple<objT, objT> checkBoxHelper(size_t l, size_t nodeIdx) {
    node* cur = &T::nodes[l][nodeIdx];

    if (T::lastLevel(l) || cur->isLeaf()) return T::getBox(l, nodeIdx);

    auto pBox = T::getBox(l, nodeIdx);
    objT pMin = std::get<0>(pBox);
    objT pMax = std::get<1>(pBox);

    for (size_t i = 0; i < 8; ++ i) {
      size_t nodeIdxNext = nodeIdx << T::dim;
      nodeIdxNext = nodeIdxNext + i;

      node* next = &T::nodes[l + 1][nodeIdxNext];

      if (next->size() > 0) {
	auto box = checkBoxHelper(l + 1, nodeIdxNext);
	objT cMin = std::get<0>(box);
	objT cMax = std::get<1>(box);

	for (int d = 0; d < 3; ++ d) {
	  EXPECT_TRUE(cMin[d] + eps>= pMin[d]);
	  EXPECT_TRUE(cMin[d] <= pMax[d] + eps);
	  EXPECT_TRUE(cMax[d] + eps>= pMin[d]);
	  EXPECT_TRUE(cMax[d] <= pMax[d] + eps);
	}
      }

    }

    return pBox;
  }

public:

  octTreeTest(parlay::slice<objT*, objT*> P, size_t _L):
    T(P, _L) {}

  void checkSum() {
    checkSumHelper(0, 0);
  }

  void checkElem() {
    checkElemHelper(0, 0);
  }

  void checkBox() {
    checkBoxHelper(0, 0);
  }
};

parlay::sequence<pargeo::point<3>> data3d(size_t n) {
  auto P = pargeo::uniformInPolyPoints<3>(n, 1);
  return std::move(P);
}

TEST(octTree_test, nodeSum) {
  using namespace pargeo;
  using pointT = pargeo::point<3>;
  auto P = data3d(10000);
  auto tree =
    octTreeTest<pointT, pointT>(make_slice(P), 7);
  tree.checkSum();
}

TEST(octTree_test, nodeElem) {
  using namespace pargeo;
  using pointT = pargeo::point<3>;
  auto P = data3d(10000);
  auto tree =
    octTreeTest<pointT, pointT>(make_slice(P), 7);
  tree.checkElem();
}

TEST(octTree_test, nodeBox) {
  using namespace pargeo;
  using pointT = pargeo::point<3>;
  auto P = data3d(10000);
  auto tree =
    octTreeTest<pointT, pointT>(make_slice(P), 7);
  tree.checkBox();
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
