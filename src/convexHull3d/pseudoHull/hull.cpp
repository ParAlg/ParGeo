#include "convexHull3d/parallelQuickHull/hull.h"
#include "convexHull3d/pseudo/hull.h"

#include <iostream>
#include <fstream>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/parlayAddon.h"

// #define PSEUDO_HULL_VERBOSE

template <class vertexT>
struct internalFacet {
  static constexpr typename vertexT::floatT numericKnob = 1e-5;

  vertexT a, b, c;

  vertexT area;

  template <class pt>
  inline typename pt::floatT signedVolume(pt d) {
    return (a-d).dot(area);
  }

  bool visible(vertexT p) {
    return signedVolume(p) > numericKnob ||
      p == a || p == b || p == c;
  }

  vertexT furthest(parlay::slice<vertexT*, vertexT*> P) {
    auto apex = vertexT();
    typename vertexT::floatT m = numericKnob;
    for (auto p: P) {
      auto m2 = signedVolume(p);
      if (m2 > m) { m = m2; apex = p; }
    }
    return apex;
  }

  vertexT furthestParallel(parlay::slice<vertexT*, vertexT*> P) {
    if (P.size() < 1000)
      return furthest(P);

    auto apex = vertexT();
    typename vertexT::floatT m = numericKnob;
    apex = parlay::max_element(P, [&](vertexT aa, vertexT bb) {
				    return signedVolume(aa) <
				      signedVolume(bb);
				  });
    return apex;
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

template<class vertexT>
parlay::sequence<vertexT>
pseudoHullHelper(internalFacet<vertexT> f,
		 parlay::sequence<vertexT> &Q) {
  using namespace parlay;
  using namespace pargeo;
  using facetT = internalFacet<vertexT>;

  if (Q.size() < 4) {
    return std::move(Q);}

  vertexT apex = f.furthestParallel(make_slice(Q)); //todo parallel

  auto interiorPt = (f.a + f.b + f.c + apex)/4;

  auto f0 = facetT(f.a-interiorPt, f.b-interiorPt, apex-interiorPt);
  auto f1 = facetT(f.b-interiorPt, f.c-interiorPt, apex-interiorPt);
  auto f2 = facetT(f.c-interiorPt, f.a-interiorPt, apex-interiorPt);
  parallel_for (0, Q.size(), [&](size_t i) {Q[i] = Q[i]-interiorPt;}, 1000);

  auto flag = sequence<int>(Q.size());
  parallel_for(0, Q.size(), [&](size_t i) {
			      if (f0.visible(Q[i])) flag[i] = 0;
			      else if (f1.visible(Q[i])) flag[i] = 1;
			      else if (f2.visible(Q[i])) flag[i] = 2;
			      else flag[i] = 3;
			    }, 1000);

  auto chunks = split_k(4, Q, flag);

  auto filtered = sequence<sequence<vertexT>>(3);

  parallel_for(0, 3, [&](size_t i) {
		       if (i == 0) {
			 filtered[0] = std::move(pseudoHullHelper<vertexT>(f0, chunks[0]));
		       } else if (i == 1) {
			 filtered[1] = std::move(pseudoHullHelper<vertexT>(f1, chunks[1]));
		       } else {
			 filtered[2] = std::move(pseudoHullHelper<vertexT>(f2, chunks[2]));
		       }
		     }, 1);

  auto Q2 = parlay::flatten(filtered);
  parallel_for (0, Q2.size(), [&](size_t i) {Q2[i] = Q2[i] + interiorPt;}, 1000);

  return std::move(Q2);
}

template<class vertexT>
parlay::sequence<vertexT>
pseudoHull(parlay::slice<vertexT*, vertexT*> P) {
  using namespace parlay;
  using namespace pargeo;
  using facetT = internalFacet<vertexT>;

  // Maximize triangle area based on fixed xMin and xMax
  size_t X[6];
  auto xx = minmax_element(P, [&](vertexT i, vertexT j) {return i[0]<j[0];});
  X[0] = xx.first - &P[0]; X[1] = xx.second - &P[0];
  auto yy = minmax_element(P, [&](vertexT i, vertexT j) {return i[1]<j[1];});
  X[2] = yy.first - &P[0]; X[3] = yy.second - &P[0];
  auto zz = minmax_element(P, [&](vertexT i, vertexT j) {return i[2]<j[2];});
  X[4] = zz.first - &P[0]; X[5] = zz.second - &P[0];

  size_t xMin, xMax;
  if (P[X[1]][0]-P[X[0]][0] > P[X[3]][1]-P[X[2]][1] && P[X[1]][0]-P[X[0]][0] > P[X[5]][2]-P[X[4]][2]) {
    xMin = X[0]; xMax = X[1];
  } else if (P[X[3]][1]-P[X[2]][1] > P[X[1]][0]-P[X[0]][0] && P[X[3]][1]-P[X[2]][1] > P[X[5]][2]-P[X[4]][2]) {
    xMin = X[2]; xMax = X[3];
  } else {
    xMin = X[4]; xMax = X[5];
  }

  vertexT x1 = P[xMin];
  vertexT x2 = P[xMax];

  auto y = max_element(P, [&](vertexT i, vertexT j) {
			    return crossProduct3d(x1-i, x2-i).length() <
			      crossProduct3d(x1-j, x2-j).length();
			  });
  size_t yApex = y - &P[0];
  vertexT y1 = P[yApex];

  // Maximize simplex volume
  vertexT area = crossProduct3d(x1-y1, x2-y1);
  auto z = max_element(P, [&](vertexT i, vertexT j) {
			    return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
			  });
  size_t zApex = z - &P[0];

  size_t c1 = xMin;
  size_t c2 = xMax;
  size_t c3 = yApex;
  size_t c4 = zApex;

  auto interiorPt = (P[c1] + P[c2] + P[c3] + P[c4])/4;

  auto Q = parlay::tabulate(P.size(), [&](size_t i) {
				       return P[i] - interiorPt;
				      });

  // Make initial facets
  auto f0 = facetT(Q[c1], Q[c2], Q[c3]);
  auto f1 = facetT(Q[c1], Q[c2], Q[c4]);
  auto f2 = facetT(Q[c3], Q[c4], Q[c2]);
  auto f3 = facetT(Q[c3], Q[c4], Q[c1]);

  auto flag = sequence<int>(Q.size());
  parallel_for(0, Q.size(), [&](size_t i) {
			      if (f0.visible(Q[i])) flag[i] = 0;
			      else if (f1.visible(Q[i])) flag[i] = 1;
			      else if (f2.visible(Q[i])) flag[i] = 2;
			      else if (f3.visible(Q[i])) flag[i] = 3;
			      else flag[i] = 4;
			    });

  auto chunks = split_k(5, Q, flag);

  auto filtered = sequence<sequence<vertexT>>(4);

  parallel_for(0, 4, [&](size_t i) {
		       if (i == 0)
			 filtered[0] = std::move(pseudoHullHelper<vertexT>(f0, chunks[0]));
		       else if (i == 1)
			 filtered[1] = std::move(pseudoHullHelper<vertexT>(f1, chunks[1]));
		       else if (i == 2)
			 filtered[2] = std::move(pseudoHullHelper<vertexT>(f2, chunks[2]));
		       else
			 filtered[3] = std::move(pseudoHullHelper<vertexT>(f3, chunks[3]));
		     }, 1);
  auto Q2 = parlay::flatten(filtered);

  return parlay::tabulate(Q2.size(), [&](size_t i) {
				       return Q2[i] + interiorPt;
				     });
}

template<class pointT>
parlay::sequence<pargeo::hull3d::facet<pointT>>
pargeo::hull3d::pseudo::compute(parlay::slice<pointT*, pointT*> P) {

#ifdef PSEUDO_HULL_VERBOSE
  timer t; t.start();
#endif

  auto Q = pseudoHull<pointT>(make_slice(P));

#ifdef PSEUDO_HULL_VERBOSE
  std::cout << "pseudohull-time = " << t.get_next() << "\n";
  std::cout << Q.size() << "\n";
#endif

  auto H = pargeo::hull3d::parallelQuickHull::compute(make_slice(Q));

#ifdef PSEUDO_HULL_VERBOSE
  std::cout << "hull-time = " << t.get_next() << "\n";
#endif

  return H;
}

template parlay::sequence<pargeo::hull3d::facet<pargeo::fpoint<3>>>
pargeo::hull3d::pseudo::compute(parlay::slice<pargeo::fpoint<3>*, pargeo::fpoint<3>*>);

template parlay::sequence<pargeo::hull3d::facet<pargeo::point<3>>>
pargeo::hull3d::pseudo::compute(parlay::slice<pargeo::point<3>*, pargeo::point<3>*>);
