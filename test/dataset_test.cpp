#include "parlay/sequence.h"
#include "pargeo/pointIO.h"
#include "dataset/uniform.h"
#include "gtest/gtest.h"

template<class Seq>
void unique(Seq& P) {
  for (size_t i = 0; i < 100; ++ i) {
    for (size_t j = i + 1; j < 100; ++ j) {
      EXPECT_TRUE(!(P[i] == P[j]));
    }
  }
}

TEST(dataset, onSphere0) {

  parlay::sequence<pargeo::point<3>> P =
    pargeo::uniformOnPolyPoints<3, pargeo::point<3>>(100, 0, 0);

  pargeo::point<3> zero;
  zero[0] = 0; zero[1] = 0; zero[2] = 0;

  for (auto p: P) EXPECT_DOUBLE_EQ(1, zero.dist(p));

  unique(P);
}

TEST(dataset, onSphereFloat0) {

  parlay::sequence<pargeo::fpoint<3>> P =
    pargeo::uniformOnPolyPoints<3, pargeo::fpoint<3>>(100, 0, 0);

  pargeo::fpoint<3> zero;
  zero[0] = 0; zero[1] = 0; zero[2] = 0;

  for (auto p: P) EXPECT_FLOAT_EQ(1, zero.dist(p));

  unique(P);
}

TEST(dataset, onSphere1) {

  parlay::sequence<pargeo::point<3>> P =
    pargeo::uniformOnPolyPoints<3, pargeo::point<3>>(100, 0, 0.1);

  pargeo::point<3> zero;
  zero[0] = 0; zero[1] = 0; zero[2] = 0;

  for (auto p: P) EXPECT_TRUE(zero.dist(p) <= 1.1 && zero.dist(p) >= 0.9);

  unique(P);
}

TEST(dataset, inSphere0) {

  parlay::sequence<pargeo::point<3>> P =
    pargeo::uniformInPolyPoints<3, pargeo::point<3>>(100, 0);

  pargeo::point<3> zero;
  zero[0] = 0; zero[1] = 0; zero[2] = 0;

  for (auto p: P) EXPECT_TRUE(zero.dist(p) <= 1.0);

  unique(P);
}

TEST(dataset, inCube0) {

  parlay::sequence<pargeo::point<3>> P =
    pargeo::uniformInPolyPoints<3, pargeo::point<3>>(100, 1);

  for (auto p: P) {
    EXPECT_TRUE(abs(p[0]) <= 1);
    EXPECT_TRUE(abs(p[1]) <= 1);
    EXPECT_TRUE(abs(p[2]) <= 1);
  }

  unique(P);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
