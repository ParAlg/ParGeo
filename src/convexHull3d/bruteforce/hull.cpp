#include <iostream>
#include <fstream>

#include "convexHull3d/bruteforce/hull.h"
#include "pargeo/getTime.h"
#include "pargeo/parlayAddon.h"
#include "pargeo/point.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"

template <typename pointT>
class internalFacet: public pargeo::hull3d::facet<pointT> {

private:

  static constexpr typename pointT::floatT eps = pointT::eps;

public:

  using baseT = pargeo::hull3d::facet<pointT>;

  pointT area;

  typename pointT::floatT signedVolume(pointT d) {
    return (baseT::a - d).dot(area);
  }

  bool visible(pointT p) {
    return signedVolume(p) > eps;
  }

  internalFacet(pointT _a, pointT _b, pointT _c):
    baseT(_a, _b, _c) {

    area = crossProduct3d(_b - _a, _c - _a);
  }

};

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::bruteforce::compute(parlay::slice<pointT*, pointT*> P) {

  parlay::sequence<pargeo::hull3d::facet<pointT>> H;

  for (size_t i = 0; i < P.size(); ++ i) {
    for (size_t j = i + 1; j < P.size(); ++ j) {
      for (size_t k = j + 1; k < P.size(); ++ k) {
	auto f = internalFacet<pointT>(P[i], P[j], P[k]);

	bool dir;
	size_t l = 0;
	for (; l < P.size(); ++ l) {
	  if (l == i || l == j || l == k) continue;
	  dir = f.visible(P[l]);
	  break;
	}

	bool onHull = true;
	for (; l < P.size(); ++ l) {
	  if (l == i || l == j || l == k) continue;
	  if (f.visible(P[l]) != dir) {
	    onHull = false;
	    break;
	  }
	}

	if (onHull) {
	  H.emplace_back(f.a, f.b, f.c);
	}

      }
    }
  }

  return H;
}

template
parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
  pargeo::hull3d::bruteforce::compute<pargeo::fpoint<3>>
  (parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*> P);

template
parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
  pargeo::hull3d::bruteforce::compute<pargeo::point<3>>
  (parlay::slice<pargeo::point<3>*, pargeo::point<3>*> P);
