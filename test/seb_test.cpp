#include "enclosingBall/ball.h"
#include "enclosingBall/sampling/seb.h"
#include "enclosingBall/welzl/seb.h"
#include "enclosingBall/scan/seb.h"

#include "parlay/primitives.h"
#include "pargeo/pointIO.h"
#include "dataset/uniform.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::point<4>> data() {
  return pargeo::uniformInPolyPoints<4, pargeo::point<4>>(1000, 0);
}

TEST(welzlMtf, compareBall) {
  auto P = data();
  pargeo::seb::ball b0 = pargeo::seb::welzl::compute(make_slice(P));
  pargeo::seb::ball b1 = pargeo::seb::welzlMtf::compute(make_slice(P));
  EXPECT_TRUE(b0 == b1);
}

TEST(welzlMtfPivot, compareBall) {
  auto P = data();
  pargeo::seb::ball b0 = pargeo::seb::welzl::compute(make_slice(P));
  pargeo::seb::ball b1 = pargeo::seb::welzlMtfPivot::compute(make_slice(P));
  EXPECT_TRUE(b0 == b1);
}

TEST(scan, compareBall) {
  auto P = data();
  pargeo::seb::ball b0 = pargeo::seb::welzl::compute(make_slice(P));
  pargeo::seb::ball b1 = pargeo::seb::scan::compute(make_slice(P));
  EXPECT_TRUE(b0 == b1);
}

TEST(sampling, compareBall) {
  auto P = data();
  pargeo::seb::ball b0 = pargeo::seb::welzl::compute(make_slice(P));
  pargeo::seb::ball b1 = pargeo::seb::sampling::compute(make_slice(P));
  EXPECT_TRUE(b0 == b1);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
