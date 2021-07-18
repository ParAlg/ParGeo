#pragma once

#include "convexHull3d/vertex.h"
//#include "convexHull3d/io.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"
#include "pargeo/parlayAddon.h"

namespace pargeo {
  namespace hull3d {

    template<class vertexT, class pointT>
    parlay::sequence<vertexT>
    point2vertexSerial(parlay::slice<pointT*, pointT*> _P) {
      parlay::sequence<vertexT> P(_P.size());
      size_t i = 0;
      for (auto _p: _P) P[i++] = vertexT(P[i].coords());
      return P;
    }

    template<class vertexT, class pointT>
    parlay::sequence<vertexT>
    point2vertex(parlay::slice<pointT*, pointT*> P) {
      return std::move(parlay::tabulate(P.size(),
				[&](size_t i) {
				  return vertexT(P[i].coords());
				}));
    }

    template<class hullT, class facetT, class pointT, class pointIn>
    hullT* initSerial(parlay::slice<pointIn*, pointIn*> PIn) {
      using vertexT = vertex<facetT, pointT>;
      using floatT = typename vertexT::floatT;

      parlay::sequence<vertexT> P = point2vertexSerial<vertexT>(PIn);

      // Maximize triangle area based on fixed xMin and xMax
      size_t X[6];
      auto xx = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[0]<j[0];});
      X[0] = xx.first - &P[0]; X[1] = xx.second - &P[0];
      auto yy = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[1]<j[1];});
      X[2] = yy.first - &P[0]; X[3] = yy.second - &P[0];
      auto zz = minmax_element_serial(P, [&](vertexT i, vertexT j) {return i[2]<j[2];});
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

      auto y = max_element_serial(P, [&](vertexT i, vertexT j) {
	  return crossProduct3d(x1-i, x2-i).length() <
	    crossProduct3d(x1-j, x2-j).length();
	});
      size_t yApex = y - &P[0];
      vertexT y1 = P[yApex];

      // Maximize simplex volume

      vertexT area = crossProduct3d(x1-y1, x2-y1);
      auto z = max_element_serial(P, [&](vertexT i, vertexT j) {
	  return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
	});
      size_t zApex = z - &P[0];

      size_t c1 = xMin;
      size_t c2 = xMax;
      size_t c3 = yApex;
      size_t c4 = zApex;

      vertexT interiorPt = (P[c1] + P[c2] + P[c3] + P[c4])/4;

      // Initialize points with visible facet link
      // auto Q = typename facetT::seqT(P.size());
      auto Q = parlay::sequence<vertexT>(P.size());

      for (size_t i=0; i<P.size(); ++i)
	Q[i] = P[i] - interiorPt; // translation

      // Make initial facets
      auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
      auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
      auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
      auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

      auto visible =
	[&](facetT *f, vertexT p) {
	  if ((f->a - p).dot(f->area) > pointT::eps)
	    return f->a != p && f->b != p && f->c != p;
	  else
	    return false;
	};

      for(size_t i=0; i<Q.size(); i++) {
	if (visible(f0, Q[i])) {
	  Q[i].seeFacet = f0;
	  f0->push_back(Q[i]);
	} else if (visible(f1, Q[i])) {
	  Q[i].seeFacet = f1;
	  f1->push_back(Q[i]);
	} else if (visible(f2, Q[i])) {
	  Q[i].seeFacet = f2;
	  f2->push_back(Q[i]);
	} else if (visible(f3, Q[i])) {
	  Q[i].seeFacet = f3;
	  f3->push_back(Q[i]);
	} else {
	  Q[i].seeFacet = nullptr;
	}
      }

      // return f0;
      hullT* linkedHull = new hullT(f0, Q, interiorPt);
      linkedHull->hSize = 4;
      linkedHull->linkFacet(f0, f1, f2, f3);
      linkedHull->linkFacet(f1, f0, f2, f3);
      linkedHull->linkFacet(f2, f1, f0, f3);
      linkedHull->linkFacet(f3, f1, f2, f0);

      return linkedHull;
    }

    template<class hullT, class facetT, class pointT, class pointIn>
    hullT* initParallel(parlay::slice<pointIn*, pointIn*> PIn) {
      using vertexT = vertex<facetT, pointT>;
      using floatT = typename vertexT::floatT;

      parlay::sequence<vertexT> P = point2vertex<vertexT>(PIn);

      // Maximize triangle area based on fixed xMin and xMax
      size_t X[6]; // extrema
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

      vertexT interiorPt = (P[c1] + P[c2] + P[c3] + P[c4])/4;

      // baseT::hSize = 4;
      // baseT::origin.setOrigin((P[c1] + P[c2] + P[c3] + P[c4])/4);

      // Initialize points with visible facet link
      // Q = sequence<vertexT>(P.size());
      auto Q = parlay::sequence<vertexT>(P.size());

      parlay::parallel_for(0, P.size(), [&](size_t i) {
				  Q[i] = P[i] - interiorPt;//translation
				});

      // Make initial facets
      auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
      auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
      auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
      auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

      auto visible =
	[&](facetT *f, vertexT p) {
	  if ((f->a - p).dot(f->area) > pointT::eps)
	    return f->a != p && f->b != p && f->c != p;
	  else
	    return false;
	};

      parlay::sequence<vertexT*> Qref = parlay::tabulate(Q.size(),
					 [&](size_t i){
					   return &Q[i];
					 });

      auto flag = parlay::sequence<int>(P.size());
      parlay::parallel_for(0, P.size(), [&](size_t i) {
				  if (visible(f0, Q[i])) {
				    flag[i] = 0; Q[i].seeFacet = f0;
				  } else if (visible(f1, Q[i])) {
				    flag[i] = 1; Q[i].seeFacet = f1;
				  } else if (visible(f2, Q[i])) {
				    flag[i] = 2; Q[i].seeFacet = f2;
				  } else if (visible(f3, Q[i])) {
				    flag[i] = 3; Q[i].seeFacet = f3;
				  } else {
				    flag[i] = 4; Q[i].seeFacet = nullptr;
				  }
				});

      auto chunks = parlay::split_k(5, Qref, flag);

      f0->reassign(chunks[0]);
      f1->reassign(chunks[1]);
      f2->reassign(chunks[2]);
      f3->reassign(chunks[3]);

      hullT* linkedHull = new hullT(f0, Q, interiorPt);
      linkedHull->hSize = 4;
      linkedHull->linkFacet(f0, f1, f2, f3);
      linkedHull->linkFacet(f1, f0, f2, f3);
      linkedHull->linkFacet(f2, f1, f0, f3);
      linkedHull->linkFacet(f3, f1, f2, f0);

      return linkedHull;
    }

  } // End namespace hull3d
} // End namespace pargeo
