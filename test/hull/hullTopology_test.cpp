#include "vertex.h"
#include "serialHull.h"
#include "convexHull3d/hull.h"
#include "convexHull3d/hullTopology.h"
#include "parlay/primitives.h"
#include "gtest/gtest.h"

TEST(hullTopology, facetTraversal) {
  sequence<vertex> P(4);
  size_t i = 0;
  P[i][0] = 0; P[i][1] = 0; P[i][2] = 0;
  i++;
  P[i][0] = 1; P[i][1] = 0; P[i][2] = 0;
  i++;
  P[i][0] = 0; P[i][1] = 1; P[i][2] = 0;
  i++;
  P[i][0] = 0; P[i][1] = 0; P[i][2] = 1;

  auto origin = pointOrigin();
  using _hull = serialHull<linkedFacet3d<vertex>, vertex, pointOrigin>;
  auto linkedHull = _hull(make_slice(P), origin);

  EXPECT_EQ(linkedHull.hullSizeDfs(), 4);

  {
  size_t count = 0;
  auto fVisit = [&](linkedFacet3d<vertex>* f) { return true; };
  auto fDo = [&](linkedFacet3d<vertex>* f) { count ++; };
  auto fStop = [&]() { return false; };
  linkedHull.dfsFacet(linkedHull.H, fVisit, fDo, fStop);

  EXPECT_EQ(count, 4);
  }

  {
  size_t count = 0;
  auto fVisit = [&](_hull::_edge e) { return true; };
  auto fDo = [&](_hull::_edge e) { count ++; };
  auto fStop = [&](){ return false; };
  linkedHull.dfsEdge(linkedHull.H, fVisit, fDo, fStop);

  EXPECT_EQ(count, 12);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
