#include "convexHull3d/bruteforce/hull.h"
#include "convexHull3d/serialQuickHull/hull.h"
#include "convexHull3d/parallelQuickHull/hull.h"
#include "convexHull3d/pseudo/hull.h"
#include "convexHull3d/sampling/hull.h"
#include "convexHull3d/divideConquer/hull.h"
#include "convexHull3d/gift/hull.h"

#include "parlay/primitives.h"
#include "pargeo/pointIO.h"
#include "dataset/uniform.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::fpoint<3>> data() {
  return pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(200, 0, 1.0);
}

bool compareFacets(parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>> H1,
		   parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>> H2) {
  using namespace pargeo;

  auto eq = [&](pargeo::fpoint<3> p1, pargeo::fpoint<3> p2) {
		return (std::abs(p1[0] - p2[0]) < p1.eps) &&
		  (std::abs(p1[1] - p2[1]) < p1.eps) &&
		  (std::abs(p1[2] - p2[2]) < p1.eps);
	      };

  auto facetEq = [&] (pargeo::hull3d::facet<pargeo::fpoint<3>> f1,
		      pargeo::hull3d::facet<pargeo::fpoint<3>> f2) {
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
  auto H1 = pargeo::hull3d::serialQuickHull::compute(make_slice(P));
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_parallel, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3d::parallelQuickHull::compute(make_slice(P));
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_divideConquer, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3d::divideConquer::compute(make_slice(P), 2);
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_pseudo, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3d::pseudo::compute(make_slice(P));
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_samping, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3d::sampling::compute(make_slice(P));
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

TEST(hull3d_gift, compareFacet) {
  auto P = data();
  auto H1 = pargeo::hull3d::gift::compute(make_slice(P));
  auto H2 = pargeo::hull3d::bruteforce::compute(make_slice(P));
  EXPECT_EQ(H1.size(), H2.size());
  EXPECT_TRUE(compareFacets(H1, H2));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
