#include "convexHull3d/bruteforce/hull.h"
#include "convexHull3d/initialize.h"
#include "convexHull3d/serialQuickHull/hull.h"
#include "quickHull.h"
#include "hullImpl.h"
#include "linkedFacet.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::serialQuickHull::compute(parlay::slice<pointT*, pointT*> P) {

  if (P.size() <= 10) return pargeo::hull3d::bruteforce::compute<pointT>(P);

  hullTopology<pointT>* linkedHull =
    pargeo::hull3d::initSerial<hullTopology<pointT>, linkedFacet<pointT>, pointT>(P);

  quickHullSerial<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>>(linkedHull);

  auto out = parlay::sequence<pargeo::hull3d::facet<pointT>>();
  linkedHull->getFacet(out);

  delete linkedHull;
  return out;
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
  pargeo::hull3d::serialQuickHull::compute<pargeo::fpoint<3>>
  (parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
  pargeo::hull3d::serialQuickHull::compute<pargeo::point<3>>
  (parlay::slice<pargeo::point<3>*, pargeo::point<3>*>);
