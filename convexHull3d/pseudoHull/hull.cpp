#include "convexHull3d/hull.h"

#include <iostream>
#include <fstream>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/parlayAddon.h"

// struct att {
//   size_t i;
// };

// using vt = pargeo::_point<3, float, float, att>;

// static std::ostream& operator<<(std::ostream& os, const vt v) {
//   for (int i=0; i<v.dim; ++i)
//     os << v.x[i] << " ";
//   return os;
// }

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
    return Q;}

  vertexT apex = f.furthest(make_slice(Q)); //todo parallel

  auto interiorPt = (f.a + f.b + f.c + apex)/4;

  /////////////////////
  auto f0 = facetT(f.a-interiorPt, f.b-interiorPt, apex-interiorPt);
  auto f1 = facetT(f.b-interiorPt, f.c-interiorPt, apex-interiorPt);
  auto f2 = facetT(f.c-interiorPt, f.a-interiorPt, apex-interiorPt);
  for (size_t i = 0; i < Q.size(); ++ i) Q[i] = Q[i]-interiorPt;

  auto flag = sequence<int>(Q.size());
  parallel_for(0, Q.size(), [&](size_t i) {
			      if (f0.visible(Q[i])) flag[i] = 0;
			      else if (f1.visible(Q[i])) flag[i] = 1;
			      else if (f2.visible(Q[i])) flag[i] = 2;
			      else flag[i] = 3;
			    });

  auto chunks = split_k_2(4, Q, flag);

  // size_t sum = 0;
  // for (auto chunk: chunks) sum += chunk.size();
  // if (sum != Q.size()) {
  //   throw std::runtime_error("sum mismatch");
  // }

  // // see if the three vertices of the facets are in the set being processed
  // bool have0 = 0;
  // bool have1 = 0;
  // bool have2 = 0;
  // bool have3 = 0;
  // for (size_t i=0; i<chunks.size()-1; ++i) {
  //   auto chunk = chunks[i];
  //   for (auto p: chunk) {
  //     if (p.attribute.i == apex.attribute.i) have0 = 1;
  //     if (p.attribute.i == f.a.attribute.i) have1 = 1;
  //     if (p.attribute.i == f.b.attribute.i) have2 = 1;
  //     if (p.attribute.i == f.c.attribute.i) have3 = 1;
  //   }
  // }

  auto filtered = sequence<sequence<vertexT>>(3);
  filtered[0] = pseudoHullHelper<vertexT>(f0, chunks[0]);
  filtered[1] = pseudoHullHelper<vertexT>(f1, chunks[1]);
  filtered[2] = pseudoHullHelper<vertexT>(f2, chunks[2]);

  auto Q2 = parlay::flatten(filtered);
  for (size_t i = 0; i < Q2.size(); ++ i)
    Q2[i] = Q2[i] + interiorPt;
  //////////////////////////////////////////////

  // bool found0 = 0;
  // bool found1 = 0;
  // bool found2 = 0;
  // bool found3 = 0;
  // for (auto p: Q2) {
  //   if (p.attribute.i == apex.attribute.i) found0 = 1;
  //   if (p.attribute.i == f.a.attribute.i) found1 = 1;
  //   if (p.attribute.i == f.b.attribute.i) found2 = 1;
  //   if (p.attribute.i == f.c.attribute.i) found3 = 1;
  // }
  // if (found0!=have0 || found1!=have1 || found2!=have2 || found3!=have3) {
  //   if (found0!=have0) std::cout << apex.attribute.i << " is missing\n";
  //   if (found1!=have1) std::cout << f.a.attribute.i << " is missing\n";
  //   if (found2!=have2) std::cout << f.b.attribute.i << " is missing\n";
  //   if (found3!=have3) std::cout << f.c.attribute.i << " is missing\n";
  //   throw std::runtime_error("facet vertices not in flattened result");
  // }

  return Q2;
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

  auto chunks = split_k_2(5, Q, flag);

  auto filtered = sequence<sequence<vertexT>>(4);
  filtered[0] = pseudoHullHelper<vertexT>(f0, chunks[0]);
  filtered[1] = pseudoHullHelper<vertexT>(f1, chunks[1]);
  filtered[2] = pseudoHullHelper<vertexT>(f2, chunks[2]);
  filtered[3] = pseudoHullHelper<vertexT>(f3, chunks[3]);

  auto Q2 = parlay::flatten(filtered);

  return parlay::tabulate(Q2.size(), [&](size_t i) {
				       return Q2[i] + interiorPt;
				     });
}

parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
pargeo::hull3dPseudo(parlay::sequence<pargeo::fpoint<3>> &P) {
  using namespace parlay;
  using pointT = pargeo::fpoint<3>;
  using floatT = pointT::floatT;
  using facetT = facet3d<pointT>;

  timer t; t.start();

  // std::ofstream myfile;
  // myfile.open("point.txt", std::ofstream::trunc);
  // for (auto p: P) {
  //   if (rand()%1000 == 0)
  //     myfile << p << "\n";
  // }
  // myfile.close();

  // sequence<vt> vtx = parlay::tabulate(P.size(), [&](size_t i){
  // 						  vt p = vt(P[i].coords());
  // 						  p.attribute.i = i;
  // 						  return p;
  // 						});
  // auto Q = pseudoHull<vt>(make_slice(vtx));

  auto Q = pseudoHull<pointT>(make_slice(P));
  std::cout << "pseudohull-time = " << t.get_next() << "\n";

  std::cout << Q.size() << "\n";

  //   std::ofstream myfilep;
  // myfilep.open("red.txt", std::ofstream::trunc);
  // for (auto p: Q) {
  //   myfilep << p << "\n";
  // }
  // myfilep.close();

  // sequence<pointT> Q2= tabulate(Q.size(), [&](size_t i){
  // 					    size_t ii = Q[i].attribute.i;
  // 					    return P[ii];
  // 					  });

  // std::cout << "filtered-pts = " << Q2.size() << "\n";

  auto H = hull3dSerial(Q);

  // auto Q3 = hull3dSerialInternal(P);
  // std::ofstream myfile;
  // myfile.open("point.txt", std::ofstream::trunc);
  // for (auto p: Q3) {
  //   bool found = false;
  //   for (auto q: Q2) {
  //     if (p[0] == q[0] &&
  // 	  p[1] == q[1] &&
  // 	  p[2] == q[2]) {
  // 	found = true;
  //     }
  //   }
  //   if (!found) myfile << p << "\n";
  // }
  // myfile.close();

  std::cout << "hull-time = " << t.get_next() << "\n";

  return H;
}
