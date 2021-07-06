#include "convexHull3d/serialQuickHull/hullImpl.h"
#include "convexHull3d/hullTopology.h"
#include "convexHull3d/vertex.h"
#include "convexHull3d/initialize.h"
#include "parlay/primitives.h"
#include "gtest/gtest.h"

TEST(hullTopology, facetTraversal) {
  using pt = pargeo::point<3>;
  using lf = pargeo::hull3d::serialQuickHull::linkedFacet<pt>;
  using vertex = pargeo::hull3d::vertex<lf, pt>;

  parlay::sequence<vertex> P(4);
  size_t i = 0;
  P[i][0] = 0; P[i][1] = 0; P[i][2] = 0;
  i++;
  P[i][0] = 1; P[i][1] = 0; P[i][2] = 0;
  i++;
  P[i][0] = 0; P[i][1] = 1; P[i][2] = 0;
  i++;
  P[i][0] = 0; P[i][1] = 0; P[i][2] = 1;

  //template<class hullT, class facetT, class pointT, class pointIn>;

  auto linkedHull = pargeo::hull3d::initSerial<
    pargeo::hull3d::serialQuickHull::hullTopology<pt>,
    lf, pt>(parlay::make_slice(P));

  EXPECT_EQ(linkedHull->hullSizeDfs(), 4);

  {
  size_t count = 0;
  auto fVisit = [&](pargeo::hull3d::serialQuickHull::linkedFacet<pt>* f) { return true; };
  auto fDo = [&](pargeo::hull3d::serialQuickHull::linkedFacet<pt>* f) { count ++; };
  auto fStop = [&]() { return false; };
  linkedHull->dfsFacet(linkedHull->H, fVisit, fDo, fStop);

  EXPECT_EQ(count, 4);
  }

  // {
  // size_t count = 0;
  // auto fVisit = [&](pargeo::hull3d::_hullTopology::_edge e) { return true; };
  // auto fDo = [&](pargeo::hull3d::_hullTopology::_edge e) { count ++; };
  // auto fStop = [&](){ return false; };
  // linkedHull->dfsEdge(linkedHull->H, fVisit, fDo, fStop);

  // EXPECT_EQ(count, 12);
  // }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
