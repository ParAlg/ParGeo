// This code is part of the Pargeo Library
// Copyright (c) 2021 Yiqiu Wang and the Pargeo Team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "convexHull3d/hullTopology.h"

#include "linkedFacet.h"

#include <set>
#include <limits>
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "parlay/primitives.h"
#include "parlay/monoid.h"
#include "parlay/hash_table.h"
#include "pargeo/algebra.h"
#include "pargeo/parlayAddon.h"

namespace pargeo {
  namespace hull3d {
    namespace parallelQuickHull {
      template <class pointT>
      class hullTopology;
    }
  }
}

template <class pointT>
class pargeo::hull3d::parallelQuickHull::hullTopology :
  public pargeo::hull3d::_hullTopology<linkedFacet<pointT>, pargeo::hull3d::vertex<linkedFacet<pointT>, pointT>> {

  using facetT = linkedFacet<pointT>;
  using vertexT = pargeo::hull3d::vertex<facetT, pointT>;
  using baseT = _hullTopology<facetT, vertexT>;

  inline bool visible(facetT* f, vertexT p) {
    return (f->a - p).dot(f->area) > baseT::numericKnob;
  }

  inline bool keep(facetT* f, vertexT p) {
    if ((f->a - p).dot(f->area) > baseT::numericKnob)
      return f->a != p && f->b != p && f->c != p;
    else
      return false;
  }

 public:
  parlay::sequence<vertexT> Q;

  hullTopology(facetT* f, parlay::sequence<vertexT>& P, vertexT interiorPt):
    baseT() {
    baseT::interiorPt = interiorPt;
    baseT::H = f;
    Q = std::move(P);
    initMap();
    updateMap();
  }

  void redistributeSerial(parlay::slice<facetT**, facetT**> facetsBeneath,
			  parlay::slice<facetT**, facetT**> newFacets) {

    baseT::hSize += newFacets.size() - facetsBeneath.size();

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numPts(); ++j) { // Point loop
	facetsBeneath[i]->pts(j)->seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (keep(newFacets[k], *facetsBeneath[i]->pts(j))) {
	    facetsBeneath[i]->pts(j)->seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->pts(j));
	    break;
	  }}}}
  }

  void redistributeParallel(parlay::slice<facetT**, facetT**> facetsBeneath,
			    parlay::slice<facetT**, facetT**> newFacets) {

    parlay::write_add(&(baseT::hSize), newFacets.size() - facetsBeneath.size());

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    auto tmpBuffer = parlay::sequence<vertexT*>(fn);
    fn = 0; // Used the second time as an offset counter
    for(int j=0; j<nf; ++j) {
      parlay::parallel_for(0, facetsBeneath[j]->numPts(),
		   [&](size_t x){tmpBuffer[fn+x] = facetsBeneath[j]->pts(x);});
      fn += facetsBeneath[j]->numPts();
    }

    auto flag = parlay::sequence<int>(fn);
    parlay::parallel_for(0, fn, [&](size_t i) {
			  flag[i] = nnf;
			  tmpBuffer[i]->seeFacet = nullptr;
			  for (int j=0; j<nnf; ++j) {
			    if (keep(newFacets[j], *tmpBuffer[i])) {
			      flag[i] = j;
			      tmpBuffer[i]->seeFacet = newFacets[j];
			      break;
			    }
			  }
			});

    auto chunks = parlay::split_k(nnf+1, tmpBuffer, flag);
    for (int j=0; j<nnf; ++j) {
      newFacets[j]->reassign(chunks[j]);
    }
  }

  vertexT* furthestApexParallel(facetT *_f=nullptr) {
    vertexT* apex = nullptr;

    auto fVisit = [&](facetT* f) {return true;};
    auto fDo = [&](facetT* f) {
		 if (f->numPts() > 0) apex = f->furthestParallel();
	       };
    auto fStop = [&]() { return !apex->isEmpty(); };
    baseT::dfsFacet(_f ? _f : baseT::H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  parlay::sequence<vertexT> furthestApexes(size_t n=-1) {
    parlay::sequence<vertexT> apexes;
    auto fVisit = [&](facetT* f) {return true;};
    auto fStop = [&]() {return apexes.size() >= n;};

    auto fFill = [&](facetT* f) {
		   if (f->numPts() > 0) {
		     auto apex = f->furthestParallel();
		     if (!apex.isEmpty()) apexes.push_back(apex);
		   }
		 };
    baseT::dfsFacet(baseT::H, fVisit, fFill, fStop);

    return apexes;
  }

  // The parallel version -- find from the vertex array
  // markers: the non-conflicting vertices do not have
  // If n not supplied, find furthest apexes of all facets
  parlay::sequence<vertexT> furthestApexesParallel(size_t n=-1) {
    using floatT = typename vertexT::floatT;
    //size_t m = parlay::num_workers() * 4; // number of random projections
    size_t m = n / 2;

    // parlay::sequence<vertexT> S(m*2);
    parlay::sequence<vertexT> S = parlay::tabulate(m*2, [&](size_t i){return vertexT();});

    size_t samples = m;
    auto projection = [&](size_t i) {
			size_t ii = parlay::hash64(i) % Q.size();
			vertexT a = Q[ii];
			vertexT b = Q[Q.size() - ii];
			vertexT ab = b - a;

			auto val = [&](vertexT p1) {
			  vertexT ap1 = p1 - a;
			  return ap1.dot(ab) / ab.dot(ab);};

			vertexT pMin = vertexT();
			vertexT pMax = vertexT();
			floatT minVal = std::numeric_limits<floatT>::max();
			floatT maxVal = std::numeric_limits<floatT>::lowest();
			// how to get to the points more efficiently?
			for (size_t j = ii; j < std::min(ii+samples*2, Q.size()); ++ j) {
			  vertexT p = Q[j];
			  if (p.seeFacet == nullptr) continue; // point not visible
			  floatT v = val(p);
			  if (v < minVal) {
			    minVal = v;
			    pMin = p;
			  }
			  if (v > maxVal) {
			    maxVal = v;
			    pMax = p;
			  }
			}

			S[i * 2] = pMin;
			S[i * 2 + 1] = pMax;
		      };

    // possibly need to deduplicate
    parlay::parallel_for(0, m, projection);// should have some parallelism here already
    auto S2 = parlay::filter(make_slice(S), [&](vertexT p){ return !p.isEmpty(); });

    // std::cout << S2.size() << "/" << S.size() << "\n";
    parlay::sort_inplace(S2, [&](vertexT const &p1, vertexT const &p2){
			       return p1.seeFacet <
				 p2.seeFacet;
			     });
    S2 = parlay::unique(S2, [&](vertexT const &p1, vertexT const &p2){
			       return p1.seeFacet ==
				 p2.seeFacet;
			     });

    if (S2.size() > 0) return std::move(S2); // make this more general
    else return std::move(furthestApexes(n));
  }

  parlay::sequence<size_t> mapping;
  size_t M;

  void initMap() {
    M = Q.size();
    mapping = parlay::tabulate(M, [&](size_t i){
			    return i;
			  });
  }

  size_t numConflict() {
    return M;
  }

  size_t updateMap() {
    auto flag = parlay::tabulate(M+1, [&](size_t i){
			   return 0;
			 });
    parlay::parallel_for(0, M, [&](size_t i) {
		   if (Q[mapping[i]].seeFacet != nullptr)
		     flag[i] = 1;
		 });

    size_t M2 = parlay::scan_inplace(make_slice(flag), parlay::addm<size_t>());
    flag[M] = M2;

    // std::cout << "total = " << M2 << "\n";

    if (M2 <= 0) return M2;

    auto map2 = parlay::sequence<size_t>(M2);
    parlay::parallel_for(0, M,
		 [&](size_t i){
		   if (flag[i] != flag[i+1]) map2[flag[i]] = mapping[i];
		 });
    mapping = std::move(map2);
    M = M2;

    // // check
    // std::cout << "check = ";
    // for (auto idx: mapping) {
    //   std::cout << idx << " ";
    //   if (Q[idx].seeFacet == nullptr) {
    // 	std::cout << "\n";
    // 	throw std::runtime_error("error!");
    //   }
    // }
    // std::cout << "checked\n";
    return M2;
  }

  // parlay::sequence<vertexT*> vertexDedup(parlay::slice<vertexT**, vertexT**> S) {
  //   parlay::sort_inplace(S, [&](vertexT const *p1, vertexT const *p2){
  // 			      return p1->seeFacet <
  // 				p2->seeFacet;
  // 			    });
  //   auto S2 = parlay::unique(S, [&](vertexT const *p1, vertexT const *p2){
  // 				  return p1->seeFacet ==
  // 				    p2->seeFacet;
  // 				});
  //   return std::move(S2);
  // }

  // // random projection
  // parlay::sequence<vertexT> furthestApexesParallel2(size_t n=-1) {
  //   using floatT = typename vertexT::floatT;

  //   size_t m = n / 4;

  //   auto S = parlay::sequence<vertexT>(m*2);

  //   auto projection = [&](size_t i) {
  // 			// todo produce more accurate projections
  // 			size_t ii = parlay::hash64(i) % Q.size();
  // 			vertexT a = Q[ii];
  // 			vertexT b = Q[Q.size() - ii];
  // 			vertexT ab = b - a;

  // 			auto val = [&](vertexT p1) {
  // 				     vertexT ap1 = p1 - a;
  // 				     return ap1.dot(ab) / ab.dot(ab);};

  // 			vertexT pMin = vertexT();
  // 			vertexT pMax = vertexT();
  // 			floatT minVal = std::numeric_limits<floatT>::max();
  // 			floatT maxVal = std::numeric_limits<floatT>::lowest();

  // 			for (size_t j = 0; j < m * 3; ++ j) {
  // 			  vertexT p = Q[mapping[(j + ii) % M]];
  // 			  floatT v = val(p);
  // 			  if (v < minVal) {
  // 			    minVal = v;
  // 			    pMin = p;
  // 			  }
  // 			  if (v > maxVal) {
  // 			    maxVal = v;
  // 			    pMax = p;
  // 			  }
  // 			}

  // 			S[i * 2] = pMin;
  // 			S[i * 2 + 1] = pMax;
  // 		      };

  //   parlay::parallel_for(0, m, projection, 1);

  //   // deduplicate
  //   auto S2 = vertexDedup(parlay::make_slice(S));

  //   //std::cout << S2.size() << "/" << S.size() << "\n";

  //   if (S2.size() > 0) return std::move(S2);
  //   else return std::move(furthestApexes(n));
  // }

  // quick hull
  parlay::sequence<vertexT*> furthestApexesParallel3(size_t n=-1) {
    using floatT = typename vertexT::floatT;

    size_t m = std::min(M, n);

    parlay::sequence<facetT*> F =
      parlay::tabulate(m, [&](size_t i){
		    return Q[mapping[i]].seeFacet;
		  });

    parlay::sort_inplace(make_slice(F),
		 [&](facetT* f1, facetT* f2){
		   return f1 < f2;
		 });


    F = std::move(parlay::unique(make_slice(F)));

    auto V = parlay::tabulate(F.size(), [&](size_t i){
				  return F[i]->furthestParallel();
				});

    // the dedup might be redundant
    return std::move(V);
    //return vertexDedup(make_slice(V));
  }

  // randomized incremental
  parlay::sequence<vertexT*> randomApexes(size_t n=-1) {
    using floatT = typename vertexT::floatT;

    size_t m = std::min(M, n);

    parlay::sequence<facetT*> F =
      parlay::tabulate(m, [&](size_t i){
		    return Q[mapping[i]].seeFacet;
		  });

    parlay::sort_inplace(make_slice(F),
		 [&](facetT* f1, facetT* f2){
		   return f1 < f2;
		 });

    F = std::move(parlay::unique(make_slice(F)));

    auto V = parlay::tabulate(F.size(), [&](size_t i){
				  return F[i]->randomApex(i);
				});

    return std::move(V);
  }

  vertexT* randomApex(facetT *_f=nullptr) {
    vertexT* apex = nullptr;

    auto fVisit = [&](facetT* f) {return true;};
    auto fDo = [&](facetT* f) {
		 if (f->numPts() > 0) apex = f->randomApex();
	       };
    auto fStop = [&]() { return !apex->isEmpty(); };
    baseT::dfsFacet(_f ? _f : baseT::H, fVisit, fDo, fStop);
    return apex;
  }

  // // If n not supplied, find furthest apexes of all facets
  // parlay::sequence<vertexT> furthestApexesWithSkip(size_t n=-1) {
  //   parlay::sequence<vertexT> apexes;
  //   auto fVisit = [&](facetT* f) {return true;};
  //   auto fStop = [&]() {return apexes.size() >= n;};

  //   parlay::hashtable<hash_pointer<facetT*>> V(baseT::hullSize(), hash_pointer<facetT*>());
  //   auto mark = [&](facetT* f) {V.insert(f);};
  //   auto visited = [&](facetT* f) {return V.find(f) != nullptr;};

  //   auto fRandPick = [&](facetT* f) {
  // 		       bool randSkip = parlay::hash64(size_t(f)) % 3;
  // 		       if (!randSkip && f->numPts() > 0) {
  // 			 mark(f);
  // 			 auto apex = f->furthestParallel();
  // 			 if (!apex.isEmpty()) {
  // 			   apexes.push_back(apex);
  // 			 }
  // 		       }
  // 		     };
  //   baseT::dfsFacet(baseT::H, fVisit, fRandPick, fStop);

  //   auto fFill = [&](facetT* f) {
  // 		   if (!visited(f) && f->numPts() > 0) {
  // 		     auto apex = f->furthestParallel();
  // 		     if (!apex.isEmpty()) {
  // 		       apexes.push_back(apex);
  // 		     }
  // 		   }
  // 		 };
  //   baseT::dfsFacet(baseT::H, fVisit, fFill, fStop);

  //   return apexes;
  // }

  facetT* facetWalk() {
    facetT* f = baseT::H;
    size_t fSize = f->numPts();

    auto fVisit = [&](facetT* _f) {return true;};
    auto fDo = [&](facetT* _f) {
		 if (_f->numPts() > fSize) {
		   fSize = _f->numPts();
		   f = _f;
		 }
	       };
    auto fStop = [&]() { return false; };
    baseT::dfsFacet(f, fVisit, fDo, fStop);
    return f;
  }

  /* Compute a frontier of edges in the clockwise order
      and facets to delete

     Meanwhile reserve the facet using min(apex.seeFacet*)
   */
  std::tuple<parlay::sequence<typename baseT::_edge>, parlay::sequence<facetT*>>
  computeFrontierAndReserve(vertexT apex) {
    using fwEdge = typename baseT::_edge;
    facetT* fVisible = apex.seeFacet;

    auto frontier = parlay::sequence<fwEdge>();
    auto facets = parlay::sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets.size(); ++i) {
			    if (f == facets.at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](fwEdge e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || visible(e.fb, apex))
		      return true;
		    else
		      return false;
		  };

    auto fDo = [&](fwEdge e) {
		 // Include the facet for deletion if visible
		 // Also reserve the facet
		 bool seeff = visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff)) {
		   facets.push_back(e.ff);
		   e.ff->reserve(apex);
		 }

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 // Also reserve invisible facet adjacent to a visible one
		 bool seefb = visible(e.fb, apex);
		 if (seefb && !seeff) {
		   frontier.emplace_back(e.a, e.b, e.ff, e.fb);
		   e.ff->reserve(apex);
		 }
	       };
    auto fStop = [&](){ return false;};

    baseT::dfsEdge(apex.seeFacet, fVisit, fDo, fStop);
    return std::make_tuple(std::move(frontier), std::move(facets));
  }

  bool confirmReservation(vertexT apex, parlay::slice<facetT**, facetT**> toDelete) {
    bool ok = true;
    for (auto f: toDelete) {
      if (!f->reserved(apex) ||
	  !f->abFacet->reserved(apex) ||
	  !f->bcFacet->reserved(apex) ||
	  !f->caFacet->reserved(apex) ) {
	ok = false;
	break;
      }
    }
    return ok;
  }

  void resetReservation(vertexT apex, parlay::slice<facetT**, facetT**> toDelete) {
    for (auto f: toDelete) {
      f->reservation = -1;
      f->abFacet->reservation = -1;
      f->bcFacet->reservation = -1;
      f->caFacet->reservation = -1;
    }
  }

  bool checkReset() {
    bool ok = true;
    auto fVisit = [&](facetT* f) { return true; };
    auto fDo = [&](facetT* f) {
		 if (f->reservation != -1) ok = false;
	       };
    auto fStop = [&](){ return false; };
    baseT::dfsFacet(baseT::H, fVisit, fDo, fStop);
    return ok;
  }

  size_t stats() {
    double numFacet = 0;
    double nonEmpty = 0;
    double totalPts = 0;
    double minPts = std::numeric_limits<double>::max();
    double maxPts = std::numeric_limits<double>::lowest();

    auto fVisit = [&](facetT* f) { return true; };

    auto fStat = [&](facetT* f) {
		   numFacet += 1;
		   totalPts += f->numPts();
		   if (f->numPts() > 0) {
		     minPts = std::min(minPts, (double)f->numPts());
		     maxPts = std::max(maxPts, (double)f->numPts());
		     nonEmpty += 1;
		   }
		 };

    auto fStop = [&](){ return false; };

    baseT::dfsFacet(baseT::H, fVisit, fStat, fStop);

    std::cout << ">>> stats\n";
    std::cout << " num-facets = " << numFacet << "\n";
    std::cout << " non-empty-facets = " << nonEmpty << "\n";
    std::cout << " total-pts = " << totalPts << "\n";
    std::cout << " min-pts = " << minPts << "\n";
    std::cout << " max-pts = " << maxPts << "\n";
    return 0;
  }

};
