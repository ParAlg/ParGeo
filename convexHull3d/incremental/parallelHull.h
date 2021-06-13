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

#include <set>
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "parlay/hash_table.h"
#include "pargeo/algebra.h"
#include "pargeo/parlayAddon.h"

template <class facetT, class vertexT, class originT>
class parallelHull : public _hullTopology<facetT, vertexT, originT> {

  using baseT = _hullTopology<facetT, vertexT, originT>;

 public:

  parallelHull(slice<vertexT*, vertexT*> P, // todo may not P in baseclass constructor
	       originT _origin):
    baseT() {
    baseT::origin = _origin;
    baseT::H = initialize(P);
  }

  facetT* initialize(slice<vertexT*, vertexT*> P) {
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

    baseT::hSize = 4;

    baseT::origin.setOrigin((P[c1] + P[c2] + P[c3] + P[c4])/4);

    // Initialize points with visible facet link
    auto Q = typename facetT::seqT(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				Q[i] = P[i] - baseT::origin.get();//translation
			      });

    // Make initial facets
    auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
    auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
    auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
    auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

    baseT::linkFacet(f0, f1, f2, f3);
    baseT::linkFacet(f1, f0, f2, f3);
    baseT::linkFacet(f2, f1, f0, f3);
    baseT::linkFacet(f3, f1, f2, f0);

    auto flag = sequence<int>(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				if (baseT::origin.keep(f0, Q[i])) {
				  flag[i] = 0; Q[i].attribute.seeFacet = f0;
				} else if (baseT::origin.keep(f1, Q[i])) {
				  flag[i] = 1; Q[i].attribute.seeFacet = f1;
				} else if (baseT::origin.keep(f2, Q[i])) {
				  flag[i] = 2; Q[i].attribute.seeFacet = f2;
				} else if (baseT::origin.keep(f3, Q[i])) {
				  flag[i] = 3; Q[i].attribute.seeFacet = f3;
				} else {
				  flag[i] = 4; Q[i].attribute.seeFacet = nullptr;
				}
			      });

    //auto chunks = split_k(5, &Q, flag);
    auto chunks = split_k_2(5, Q, flag);

    f0->reassign(chunks[0]);
    f1->reassign(chunks[1]);
    f2->reassign(chunks[2]);
    f3->reassign(chunks[3]);

    return f0;
  }

  void redistributeSerial(slice<facetT**, facetT**> facetsBeneath,
			  slice<facetT**, facetT**> newFacets) {

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
	facetsBeneath[i]->pts(j).attribute.seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (baseT::origin.keep(newFacets[k], facetsBeneath[i]->pts(j))) {
	    facetsBeneath[i]->pts(j).attribute.seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->pts(j));
	    break;
	  }}}}
  }

  void redistributeParallel(slice<facetT**, facetT**> facetsBeneath,
			    slice<facetT**, facetT**> newFacets) {

    parlay::write_add(&(baseT::hSize), newFacets.size() - facetsBeneath.size());

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }

    auto tmpBuffer = typename facetT::seqT(fn);
    fn = 0; // Used the second time as an offset counter
    for(int j=0; j<nf; ++j) {
      parallel_for(0, facetsBeneath[j]->numPts(),
		   [&](size_t x){tmpBuffer[fn+x] = facetsBeneath[j]->pts(x);});
      fn += facetsBeneath[j]->numPts();
    }

    auto flag = sequence<int>(fn);
    parallel_for(0, fn, [&](size_t i) {
			  flag[i] = nnf;
			  tmpBuffer[i].attribute.seeFacet = nullptr;
			  for (int j=0; j<nnf; ++j) {
			    if (baseT::origin.keep(newFacets[j], tmpBuffer[i])) {
			      flag[i] = j;
			      tmpBuffer[i].attribute.seeFacet = newFacets[j];
			      break;
			    }
			  }
			});

    //auto chunks = split_k(nnf+1, &tmpBuffer, flag);
    auto chunks = split_k_2(nnf+1, tmpBuffer, flag);
    for (int j=0; j<nnf; ++j) {
      newFacets[j]->reassign(chunks[j]);
    }
  }

  vertexT furthestApexParallel(facetT *_f=nullptr) {
    vertexT apex = vertexT();

    auto fVisit = [&](facetT* f) {return true;};
    auto fDo = [&](facetT* f) {
		 if (f->numPts() > 0) apex = f->furthestParallel();
	       };
    auto fStop = [&]() { return !apex.isEmpty(); };
    baseT::dfsFacet(_f ? _f : baseT::H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertexT> furthestApexes(size_t n=-1) {
    sequence<vertexT> apexes;
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

  // If n not supplied, find furthest apexes of all facets
  sequence<vertexT> randomApexes(size_t n=-1) {
    sequence<vertexT> apexes;
    auto fVisit = [&](facetT* f) {return true;};
    auto fStop = [&]() {return apexes.size() >= n;};

    auto fFill = [&](facetT* f) {
		   if (f->numPts() > 0) {
		     auto apex = f->pts(parlay::hash64(f->numPts()) % f->numPts());
		     if (!apex.isEmpty()) apexes.push_back(apex);
		   }
		 };
    baseT::dfsFacet(baseT::H, fVisit, fFill, fStop);

    return apexes;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertexT> furthestApexesWithSkip(size_t n=-1) {
    sequence<vertexT> apexes;
    auto fVisit = [&](facetT* f) {return true;};
    auto fStop = [&]() {return apexes.size() >= n;};

    hashtable<hash_pointer<facetT*>> V(baseT::hullSize(), hash_pointer<facetT*>());
    auto mark = [&](facetT* f) {V.insert(f);};
    auto visited = [&](facetT* f) {return V.find(f) != nullptr;};

    auto fRandPick = [&](facetT* f) {
		       bool randSkip = parlay::hash64(size_t(f)) % 3;
		       if (!randSkip && f->numPts() > 0) {
			 mark(f);
			 auto apex = f->furthestParallel();
			 if (!apex.isEmpty()) {
			   apexes.push_back(apex);
			 }
		       }
		     };
    baseT::dfsFacet(baseT::H, fVisit, fRandPick, fStop);

    auto fFill = [&](facetT* f) {
		   if (!visited(f) && f->numPts() > 0) {
		     auto apex = f->furthestParallel();
		     if (!apex.isEmpty()) {
		       apexes.push_back(apex);
		     }
		   }
		 };
    baseT::dfsFacet(baseT::H, fVisit, fFill, fStop);

    return apexes;
  }

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

     Meanwhile reserve the facet using min(apex.attribute.seeFacet*)
   */
  //tuple<sequence<typename baseT::_edge>*, sequence<facetT*>*> computeFrontierAndReserve(vertexT apex) {
  tuple<sequence<typename baseT::_edge>, sequence<facetT*>> computeFrontierAndReserve(vertexT apex) {
    using fwEdge = typename baseT::_edge;
    facetT* fVisible = apex.attribute.seeFacet;

    // auto frontier = new sequence<fwEdge>();
    // auto facets = new sequence<facetT*>();
    auto frontier = sequence<fwEdge>();
    auto facets = sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets.size(); ++i) {
			    if (f == facets.at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](fwEdge e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || baseT::origin.visible(e.fb, apex))
		      return true;
		    else
		      return false;
		  };

    auto fDo = [&](fwEdge e) {
		 // Include the facet for deletion if visible
		 // Also reserve the facet
		 bool seeff = baseT::origin.visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff)) {
		   facets.push_back(e.ff);
		   e.ff->reserve(apex);
		 }

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 // Also reserve invisible facet adjacent to a visible one
		 bool seefb = baseT::origin.visible(e.fb, apex);
		 if (seefb && !seeff) {
		   frontier.emplace_back(e.a, e.b, e.ff, e.fb);
		   e.ff->reserve(apex);
		 }
	       };
    auto fStop = [&](){ return false;};

    baseT::dfsEdge(apex.attribute.seeFacet, fVisit, fDo, fStop);
    return make_tuple(frontier, facets);
  }

  bool confirmReservation(vertexT apex, slice<facetT**, facetT**> toDelete) {
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

  void resetReservation(vertexT apex, slice<facetT**, facetT**> toDelete) {
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
		     minPts = min(minPts, (double)f->numPts());
		     maxPts = max(maxPts, (double)f->numPts());
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
