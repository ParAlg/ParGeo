#include "convexHull3d/hull.h"
#include "convexHull3d/bruteforceHull.h"
#include "convexHull3d/serialHull.h"
#include "convexHull3d/pseudoHull.h"
#include "convexHull3d/gridHull.h"
#include "convexHull3d/samplingHull.h"
#include "convexHull3d/searchHull.h"
#include "convexHull3d/concurrentHull.h"
#include "convexHull3d/incrementalHull.h"

#include "parlay/primitives.h"
#include "pargeo/pointIO.h"
#include "dataset/uniform.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::fpoint<3>> data() {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(100, 0);
}

bool compareFacets(parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>> H1,
		   parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>> H2) {
  using namespace pargeo;

  double knob = 1e-6;

  auto eq = [&](pargeo::fpoint<3> p1, pargeo::fpoint<3> p2) {
		return (std::abs(p1[0] - p2[0]) < knob) &&
		  (std::abs(p1[1] - p2[1]) < knob) &&
		  (std::abs(p1[2] - p2[2]) < knob);
	      };

  auto facetEq = [&] (facet3d<pargeo::fpoint<3>> f1,
		      facet3d<pargeo::fpoint<3>> f2) {
		   return
		     (eq(f1.a,f2.a) && eq(f1.b,f2.b) && eq(f1.c,f2.c)) ||
		     (eq(f1.a,f2.a) && eq(f1.b,f2.c) && eq(f1.c,f2.b)) ||
		     (eq(f1.a,f2.b) && eq(f1.b,f2.a) && eq(f1.c,f2.c)) ||
		     (eq(f1.a,f2.b) && eq(f1.b,f2.c) && eq(f1.c,f2.a)) ||
		     (eq(f1.a,f2.c) && eq(f1.b,f2.a) && eq(f1.c,f2.b)) ||
		     (eq(f1.a,f2.c) && eq(f1.b,f2.b) && eq(f1.c,f2.a));
		 };

  for(auto f1: H1) {
    bool found = false;
    for(auto f2: H2) {
      if (facetEq(f1, f2)) found = true;
    }
    if (found) continue;
    else {
      for(auto f: H1) std::cout << f.a << ", " << f.b << ", " << f.c << "\n";
      std::cout << "---\n";
      for(auto f: H2) std::cout << f.a << ", " << f.b << ", " << f.c << "\n";
      std::cout << "Facet not found " << f1.a << ", " << f1.b << ", " << f1.c << "\n";
      return false;
    }
  }
  return true;
}

TEST(hull3d_serial, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dSerial(P);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_incremental, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dIncremental(P, 2);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_concurrent, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dConcurrent(P, 2);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_pseudo, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dPseudo(P);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_grid, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dGrid(P, 2, false);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_samping, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dSampling(P, 0.01);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_search, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3dSearch(P);
  auto H2 = pargeo::hull3dBruteforce(P);
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
