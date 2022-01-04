#include "convexHull3d/bruteforce/hull.h"

#include <iostream>
#include <fstream>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/parlayAddon.h"

template <class vertexT>
struct internalFacet {
  static constexpr typename vertexT::floatT numericKnob = vertexT::eps;

  vertexT a, b, c;

  vertexT area;

  template <class pt>
  inline typename pt::floatT signedVolume(pt d) {
    return (a-d).dot(area);
  }

  bool visible(vertexT p) {
    return signedVolume(p) > numericKnob;
  }

  internalFacet(vertexT _a, vertexT _b, vertexT _c):
    a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > numericKnob)
      std::swap(b, c);
    area = crossProduct3d(b-a, c-a);
  }

  ~internalFacet() {
  }

};

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::bruteforce::compute(parlay::slice<pointT*, pointT*> P) {
  using namespace parlay;
  using floatT = typename pointT::floatT;
  using facetT = pargeo::hull3d::facet<pointT>;
  using fc = internalFacet<pointT>;

  sequence<facetT> H;

  for (size_t i = 0; i < P.size(); ++ i) {
    for (size_t j = i + 1; j < P.size(); ++ j) {
      for (size_t k = j + 1; k < P.size(); ++ k) {
	auto f = fc(P[i], P[j], P[k]);

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
