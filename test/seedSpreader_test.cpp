#include "parlay/sequence.h"
#include "dataset/seedSpreader.h"
#include "gtest/gtest.h"
#include "pargeo/pointIO.h"
#include "pargeo/point.h"

TEST(dataset, simden0) {
  static const int dim = 2;
  auto P = pargeo::seedSpreader::simdenGenerator<dim>(1000, 0.01);
  // pargeo::pointIO::writePointsToFile(P, "simden.txt");
  EXPECT_EQ(P.size(), 1000);
}

TEST(dataset, varden0) {
  static const int dim = 2;
  auto P = pargeo::seedSpreader::vardenGenerator<dim>(1000, 0.01);
  // pargeo::pointIO::writePointsToFile(P, "varden.txt");
  EXPECT_EQ(P.size(), 1000);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
