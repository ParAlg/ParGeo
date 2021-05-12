#include "pargeo/pointIO.h"
#include "pargeo/closestPair.h"
#include "gtest/gtest.h"

TEST(closestPair_divideConquer, test2d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "../datasets/2d_100.txt";
  int dim = readHeader(filePath);
  parlay::sequence<pargeo::point<2>> P =
    readPointsFromFile<pargeo::point<2>>(filePath);

  auto cp1 = closestPairDC<2>(P);
  auto cp2 = bruteForceParallel<2>(make_slice(P));
  if (cp1.u != cp2.u) swap(cp1.u, cp1.v);
  EXPECT_EQ(cp1.dist, cp2.dist);
  EXPECT_EQ(cp1.u[0], cp2.u[0]);
  EXPECT_EQ(cp1.u[1], cp2.u[1]);
  EXPECT_EQ(cp1.v[0], cp2.v[0]);
  EXPECT_EQ(cp1.v[1], cp2.v[1]);
}

TEST(closestPair_divideConquer, test3d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "../datasets/3d_40.txt";
  int dim = readHeader(filePath);
  parlay::sequence<pargeo::point<3>> P =
    readPointsFromFile<pargeo::point<3>>(filePath);

  auto cp1 = closestPairDC<3>(P);
  auto cp2 = bruteForceParallel<3>(make_slice(P));
  if (cp1.u != cp2.u) swap(cp1.u, cp1.v);
  EXPECT_EQ(cp1.dist, cp2.dist);
  EXPECT_EQ(cp1.u[0], cp2.u[0]);
  EXPECT_EQ(cp1.u[1], cp2.u[1]);
  EXPECT_EQ(cp1.u[2], cp2.u[2]);
  EXPECT_EQ(cp1.v[0], cp2.v[0]);
  EXPECT_EQ(cp1.v[1], cp2.v[1]);
  EXPECT_EQ(cp1.v[2], cp2.v[2]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
