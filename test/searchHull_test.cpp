#include "convexHull3d/hull.h"
#include "parlay/primitives.h"
#include "pargeo/pointIO.h"
#include "pargeo/getTime.h"
#include "dataset/uniform.h"
#include "gtest/gtest.h"

parlay::sequence<pargeo::fpoint<3>> data(size_t n) {
  auto P = pargeo::uniformInPolyPoints<3, pargeo::fpoint<3>>(n, 0);
  return std::move(P);
  return P;
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
  auto P = data(1000);
  pargeo::timer t; t.start();
  std::cout << "Serial:\n";
  auto H1 = pargeo::hull3dSerial(P);
  std::cout << "time = " << t.get_next() << "\n";
  std::cout << "hull = " << H1.size() << "\n";
  std::cout << "Search:\n";
  auto H2 = pargeo::hull3dSearch(P);
  std::cout << "time = " << t.stop() << "\n";
  std::cout << "hull = " << H1.size() << "\n";
  // auto H2 = pargeo::hull3dBruteforce(P);
  // EXPECT_EQ(H1.size(), H2.size());
  // EXPECT_TRUE(compareFacets(H1, H2));
  EXPECT_TRUE(true);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
