#include "pargeo/pointIO.h"
#include "closestPair/closestPair.h"
#include "gtest/gtest.h"

using namespace std;
using namespace parlay;
using namespace pargeo;

template<int dim>
pair<point<dim>, point<dim>> bruteForce(slice<point<dim>*,point<dim>*> P) {
  using floatT = typename point<dim>::floatT;

  point<dim> u = P[0];
  point<dim> v = P[1];
  floatT dist = u.dist(v);

  for(size_t i = 0; i < P.size(); ++ i) {
    for(size_t j = i+1; j < P.size(); ++ j) {
      floatT distt = P[i].dist(P[j]);
      if (distt < dist) {
	u = P[i];
	v = P[j];
	dist = distt;
      }
    }
  }

  return {u, v};
}

TEST(closestPair_divideConquer, test2d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/2d_100.txt";
  int dim = readHeader(filePath);
  parlay::sequence<pargeo::point<2>> P =
    readPointsFromFile<pargeo::point<2>>(filePath);

  auto cp1 = closestPair<2>(P);
  auto cp2 = bruteForce<2>(make_slice(P));
  if (cp1.first != cp2.first) swap(cp1.first, cp1.second);
  EXPECT_EQ(cp1.first[0], cp2.first[0]);
  EXPECT_EQ(cp1.first[1], cp2.first[1]);
  EXPECT_EQ(cp1.second[0], cp2.second[0]);
  EXPECT_EQ(cp1.second[1], cp2.second[1]);
}

TEST(closestPair_divideConquer, test3d) {
  using namespace pargeo;
  using namespace pargeo::pointIO;

  auto filePath = "datasets/3d_40.txt";
  int dim = readHeader(filePath);
  parlay::sequence<pargeo::point<3>> P =
    readPointsFromFile<pargeo::point<3>>(filePath);

  auto cp1 = closestPair<3>(P);
  auto cp2 = bruteForce<3>(make_slice(P));
  if (cp1.first != cp2.first) swap(cp1.first, cp1.second);
  EXPECT_EQ(cp1.first[0], cp2.first[0]);
  EXPECT_EQ(cp1.first[1], cp2.first[1]);
  EXPECT_EQ(cp1.first[2], cp2.first[2]);
  EXPECT_EQ(cp1.second[0], cp2.second[0]);
  EXPECT_EQ(cp1.second[1], cp2.second[1]);
  EXPECT_EQ(cp1.second[2], cp2.second[2]);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
