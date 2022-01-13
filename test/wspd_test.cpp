#include "pargeo/pointIO.h"
#include "kdTree/kdTree.h"
#include "gtest/gtest.h"

template <int dim>
inline void testWspd(parlay::sequence<pargeo::point<dim>>& P,
		     parlay::sequence<pargeo::kdTree::wsp<pargeo::kdTree::node<dim, pargeo::point<dim>>>> pairs,
		     double constant) {
  using namespace std;
  using namespace pargeo;
  using pointT = point<dim>;
  using nodeT = kdTree::node<dim, point<dim>>;

  // Check that each pair is well separated
  for (auto pair: pairs) {
    EXPECT_TRUE(geomWellSeparated(pair.u, pair.v, constant));
  }

  // Check that each pair of point is in a pair
  auto find = [&](pointT p1, pointT p2){
		for (auto pair: pairs) {

		  for (size_t i = 0; i < pair.u->size(); ++ i) {
		    for (size_t j = 0; j < pair.v->size(); ++ j) {

		      if (p1 == *pair.u->at(i) && p2 == *pair.v->at(j))
			return true;
		      if (p2 == *pair.u->at(i) && p1 == *pair.v->at(j))
			return true;
		    }
		  }

		}
		return false;
		abort();
	      };

  for (size_t i = 0; i < P.size(); ++ i) {
    for (size_t j = i + 1; j < P.size(); ++ j) {
      EXPECT_TRUE(find(P[i], P[j]));
    }
  }

}

TEST(wspd_test, testSerial2d) {
  using namespace parlay;
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/2d_100.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<2>> P =
    readPointsFromFile<pargeo::point<2>>(filePath);

  kdTree::node<2, point<2>>* tree = kdTree::build<2, point<2>>(P, true, 1);

  auto pairs1 = kdTree::wspdSerial(tree, 2);

  testWspd<2>(P, pairs1, 2);

  auto pairs2 = kdTree::wspdSerial(tree, 2.2);

  testWspd<2>(P, pairs2, 2.2);
}

TEST(wspd_test, testParallel2d) {
  using namespace parlay;
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/2d_100.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<2>> P =
    readPointsFromFile<pargeo::point<2>>(filePath);

  kdTree::node<2, point<2>>* tree = kdTree::build<2, point<2>>(P, true, 1);

  auto pairs1 = kdTree::wellSeparatedPairDecomp(tree, 2);

  testWspd<2>(P, pairs1, 2);

  auto pairs2 = kdTree::wellSeparatedPairDecomp(tree, 3.9);

  testWspd<2>(P, pairs2, 3.9);
}

TEST(wspd_test, testSerial3d) {
  using namespace parlay;
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/3d_40.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<3>> P =
    readPointsFromFile<pargeo::point<3>>(filePath);

  kdTree::node<3, point<3>>* tree = kdTree::build<3, point<3>>(P, true, 1);

  auto pairs1 = kdTree::wspdSerial(tree, 2);

  testWspd<3>(P, pairs1, 2);

  auto pairs2 = kdTree::wspdSerial(tree, 2.2);

  testWspd<3>(P, pairs2, 2.2);
}

TEST(wspd_test, testParallel3d) {
  using namespace parlay;
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/3d_40.txt";
  int dim = readHeader(filePath);

  parlay::sequence<pargeo::point<3>> P =
    readPointsFromFile<pargeo::point<3>>(filePath);

  kdTree::node<3, point<3>>* tree = kdTree::build<3, point<3>>(P, true, 1);

  auto pairs1 = kdTree::wellSeparatedPairDecomp(tree, 2);

  testWspd<3>(P, pairs1, 2);

  auto pairs2 = kdTree::wellSeparatedPairDecomp(tree, 3.9);

  testWspd<3>(P, pairs2, 3.9);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
