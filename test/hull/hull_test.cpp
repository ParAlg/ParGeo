#include "convexHull3d/hull.h"
#include "pargeo/pointIO.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::fpoint<3>> data40() {
  auto filePath = "datasets/3d_40.txt";
  int dim = pargeo::pointIO::readHeader(filePath);

  parlay::sequence<pargeo::fpoint<3>> P =
    pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);

  return P;
}

TEST(hull3d_serial, hullSize) {
  auto P = data40();
  auto H1 = pargeo::hull3dSerial(P);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
}

TEST(hull3d_incremental, hullSize) {
  auto P = data40();
  auto H1 = pargeo::hull3dIncremental(P, 2);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
}

TEST(hull3d_concurrent, hullSize) {
  auto P = data40();
  auto H1 = pargeo::hull3dConcurrent(P, 2);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
}

TEST(hull3d_pseudo, hullSize) {
  auto P = data40();
  auto H1 = pargeo::hull3dPseudo(P);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
}

TEST(hull3d_grid, hullSize) {
  auto P = data40();
  auto H1 = pargeo::hull3dGrid(P, 2, false);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
