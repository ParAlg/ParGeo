// This code is part of the Pargeo Library
// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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

#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

#include <atomic>
#include <vector>
#include <stack>
#include <tuple>
#include "parlay/sequence.h"
#include "parlay/hash_table.h"
#include "geometry/point.h"
#include "geometry/algebra.h"
#include "common/get_time.h"
#include "split.h"

using namespace std;
using namespace parlay;
using namespace parlay::internal;

// Example for hashing numeric values.
// T must be some integer type
template <class T>
struct hash_pointer {
  using eType = T;
  using kType = T;
  eType empty() { return nullptr; }
  kType getKey(eType v) { return v; }
  size_t hash(kType v) { return static_cast<size_t>(hash64(size_t(v))); }
  int cmp(kType v, kType b) { return (v > b) ? 1 : ((v == b) ? 0 : -1); }
  bool replaceQ(eType, eType) { return 0; }
  eType update(eType v, eType) { return v; }
  bool cas(eType* p, eType o, eType n) {
    return std::atomic_compare_exchange_strong_explicit(
      reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
  }
};

template <class facetT, class vertexT, class originT>
struct _hull {

private:

  // A linked structure for the facets
  facetT* H;

  // The number of facets in H
  std::atomic<size_t> hSize;

  /* Depth-first hull traversal (no facet repeat)
  */
  template <class F, class G, class H>
  void dfsFacet(facetT* start, F& fVisit, G& fDo, H& fStop) {
    
    hashtable<hash_pointer<facetT*>> V(hSize, hash_pointer<facetT*>());
    auto mark = [&](facetT* f) {V.insert(f);};
    auto visited = [&](facetT* f) {return V.find(f) != nullptr;};

    stack<_edge> S;

    S.emplace(start->b, start->a, start, nullptr);

    while (S.size() > 0) {
      _edge e = S.top(); S.pop();
      if (!visited(e.ff) && fVisit(e)) {
	fDo(e);
	mark(e.ff);

	S.emplace(e.ff->a, e.ff->b, e.ff->abFacet, e.ff);
	S.emplace(e.ff->c, e.ff->a, e.ff->caFacet, e.ff);
	S.emplace(e.ff->b, e.ff->c, e.ff->bcFacet, e.ff);
      }
      if (fStop()) break;
    }
  }

  /* Clockwise, depth-first hull traversal, visits each oriented edge once
     - (facetT*) start: starting facet
     - (func _edge -> bool) fVisit : whether to visit the target facet of an advancing edge
     - (func _edge -> void) fDo : what do to with the facet in the advancing edge
     - (func void -> bool)  fStop : whether to stop
  */
  template <class F, class G, class H>
  void dfsEdge(facetT* start, F& fVisit, G& fDo, H& fStop) {
    using edgeT = _edge;

    // Quadratic iteration of V seems to be fast
    // as the involved facets are few
    sequence<edgeT> V;
    auto mark = [&](edgeT f) {V.push_back(f);};
    auto visited = [&](edgeT f) {
		     for (auto g: V) {
		       if (f == g) return true;
		     }
		     return false;};

    stack<edgeT> S;

    /* Create initial advancing edge (b,a)
         o start->c
        / \ ---> start
       o---o
start->b  start->a
    */
    S.emplace(start->b, start->a, start, nullptr);
    while (S.size() > 0) {
      _edge e = S.top(); S.pop();
      if (!visited(e) && fVisit(e)) {
	fDo(e);
	mark(e);
	/* Push in ccw order, pop in cw order; start from (e.ff.a, e.ff.b)
	   e.ff.b
	     o
	    / \ ---> e.ff
	   o---o
e.b==e.ff.a    e.a==e.ff.c
	*/
	if (e.ff->a == e.b) {
	  S.emplace(e.ff->a, e.ff->b, e.ff->abFacet, e.ff);
	  S.emplace(e.ff->c, e.ff->a, e.ff->caFacet, e.ff);
	  S.emplace(e.ff->b, e.ff->c, e.ff->bcFacet, e.ff);
	} else if (e.ff->b == e.b) {
	  S.emplace(e.ff->b, e.ff->c, e.ff->bcFacet, e.ff);
	  S.emplace(e.ff->a, e.ff->b, e.ff->abFacet, e.ff);
	  S.emplace(e.ff->c, e.ff->a, e.ff->caFacet, e.ff);
	} else if (e.ff->c == e.b) {
	  S.emplace(e.ff->c, e.ff->a, e.ff->caFacet, e.ff);
	  S.emplace(e.ff->b, e.ff->c, e.ff->bcFacet, e.ff);
	  S.emplace(e.ff->a, e.ff->b, e.ff->abFacet, e.ff);
	}
      }
      if (fStop()) break;
    }
  }

public:
  // An arbitrary coordinate class located within the hull
  // it also contains primitives for visibility test
  originT origin;

  // Link f with ab, bc, ca; the edge matching is automatic -- input in any order
  void linkFacet(facetT* f, facetT* ab, facetT* bc, facetT* ca) {
    using fc = facetT;
    fc* F[3]; F[0]=ab; F[1]=bc; F[2]=ca;
    auto findFacet = [&](vertexT v1, vertexT v2) {
		       for(int i=0; i<3; ++i) {
			 if ((F[i]->a==v1 && F[i]->b==v2) || (F[i]->b==v1 && F[i]->a==v2) ||
			     (F[i]->b==v1 && F[i]->c==v2) || (F[i]->c==v1 && F[i]->b==v2) ||
			     (F[i]->c==v1 && F[i]->a==v2) || (F[i]->a==v1 && F[i]->c==v2)) {
			   return F[i];
			 }
		       }
#ifdef VERBOSE
		       cout << "Facets: " << endl;
		       cout << "linking " << *f << " with " << endl;
		       cout << *ab << endl;
		       cout << *bc << endl;
		       cout << *ca << endl;
		       for(int i=0; i<3; ++i)
			 cout << F[i]->a.attribute.i << " "
			      << F[i]->b.attribute.i << " "
			      << F[i]->c.attribute.i << endl;
		       cout << "match v1 = " << v1.attribute.i << endl;
		       cout << "match v2 = " << v2.attribute.i << endl;
#endif
		       throw std::runtime_error("Facet linking failure.");
		     };

    auto linkEdge = [&](fc* f1, fc* f2, vertexT v1, vertexT v2) {
		      if ( (f2->a==v1 && f2->b==v2) || (f2->a==v2 && f2->b==v1) ) {
			f2->abFacet = f1;
		      } else if ( (f2->b==v1 && f2->c==v2) || (f2->b==v2 && f2->c==v1) ) {
			f2->bcFacet = f1;
		      } else if ( (f2->c==v1 && f2->a==v2) || (f2->c==v2 && f2->a==v1) ) {
			f2->caFacet = f1;
		      } else {
			throw std::runtime_error("Edge linking failure.");
		      }
		    };

    f->abFacet = findFacet(f->a, f->b);
    linkEdge(f, f->abFacet, f->a, f->b);
    f->bcFacet = findFacet(f->b, f->c);
    linkEdge(f, f->bcFacet, f->b, f->c);
    f->caFacet = findFacet(f->c, f->a);
    linkEdge(f, f->caFacet, f->c, f->a);
  }
  
  /* An edge sub structure for advancing direction in a traversal

     ff: forward facet to advance to
     fb: parent facet from which we come

        o
       / \ --->ff
    a o---o b
       \ / --->fb
        o
  */
  struct _edge {
    vertexT a, b;
    facetT* ff;
    facetT* fb;

    _edge(vertexT _a, vertexT _b, facetT* _ff, facetT* _fb) {
      a = _a; b = _b; ff = _ff; fb = _fb;}

    friend bool operator==(_edge e1, _edge e2) {
      return e1.a==e2.a && e1.b==e2.b;}//todo
  };

  void setStart(facetT* _H) {H = _H;}

  _hull(slice<vertexT*, vertexT*> P, originT _origin) {
    origin = _origin;

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
  
    hSize = 4;

    origin.setOrigin((P[c1] + P[c2] + P[c3] + P[c4])/4);

    // Initialize points with visible facet link
    auto Q = typename facetT::seqT(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				Q[i] = P[i] - origin.get();//translation
			      });

    // Make initial facets
    auto f0 = new facetT(Q[c1], Q[c2], Q[c3]);
    auto f1 = new facetT(Q[c1], Q[c2], Q[c4]);
    auto f2 = new facetT(Q[c3], Q[c4], Q[c2]);
    auto f3 = new facetT(Q[c3], Q[c4], Q[c1]);

    linkFacet(f0, f1, f2, f3);
    linkFacet(f1, f0, f2, f3);
    linkFacet(f2, f1, f0, f3);
    linkFacet(f3, f1, f2, f0);

#ifdef SERIAL
    bool serial = true;
#else
    bool serial = false;
#endif

    if (Q.size() < 1000 || serial) {

      for(size_t i=0; i<P.size(); i++) {
	if (origin.keep(f0, Q[i])) {
	  Q[i].attribute.seeFacet = f0;
	  f0->push_back(Q[i], &origin);
	} else if (origin.keep(f1, Q[i])) {
	  Q[i].attribute.seeFacet = f1;
	  f1->push_back(Q[i], &origin);
	} else if (origin.keep(f2, Q[i])) {
	  Q[i].attribute.seeFacet = f2;
	  f2->push_back(Q[i], &origin);
	} else if (origin.keep(f3, Q[i])) {
	  Q[i].attribute.seeFacet = f3;
	  f3->push_back(Q[i], &origin);
	} else {
	  Q[i].attribute.seeFacet = nullptr;
	}
      }

    } else {

      auto flag = sequence<int>(P.size());
      parallel_for(0, P.size(), [&](size_t i) {
				  if (origin.keep(f0, Q[i])) {
				    flag[i] = 0; Q[i].attribute.seeFacet = f0;
				  } else if (origin.keep(f1, Q[i])) {
				    flag[i] = 1; Q[i].attribute.seeFacet = f1;
				  } else if (origin.keep(f2, Q[i])) {
				    flag[i] = 2; Q[i].attribute.seeFacet = f2;
				  } else if (origin.keep(f3, Q[i])) {
				    flag[i] = 3; Q[i].attribute.seeFacet = f3;
				  } else {
				    flag[i] = 4; Q[i].attribute.seeFacet = nullptr;
				  }
				});

      auto chunks = split_k(5, &Q, flag);

      f0->reassign(chunks[0], &origin);
      f1->reassign(chunks[1], &origin);
      f2->reassign(chunks[2], &origin);
      f3->reassign(chunks[3], &origin);
    }

#ifdef VERBOSE
    cout << "initial-hull:" << endl;
    cout << " facet 0, area = " << f0->area.length()/2 << " #pts = " << f0->numPts() << endl;
    cout << " facet 1, area = " << f1->area.length()/2 << " #pts = " << f1->numPts() << endl;
    cout << " facet 2, area = " << f2->area.length()/2 << " #pts = " << f2->numPts() << endl;
    cout << " facet 3, area = " << f3->area.length()/2 << " #pts = " << f3->numPts() << endl;
    cout << " volume = " << abs((f0->a-Q[c4]).dot(f0->area)/6) << endl;
#endif

    H = f0;
  }

  /* Choose a random outside vertex visible to some facet
   */
  vertexT randomApex() {
    if (H->size() > 0) return H->at(0);
    vertexT apex = vertexT();
    auto fVisit = [&](_edge e) {return true;};
    auto fDo = [&](_edge e) {
		 if (e.ff->numVisiblePts() > 0)
		   apex = e.ff->visiblePts(0);};
    auto fStop = [&]() {
		   if (!apex.isEmpty()) return true;
		   else return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertexT> randomApexes(size_t n=-1) {
    sequence<vertexT> apexes;

    auto fVisit = [&](_edge e) {return true;};
    auto fDo = [&](_edge e) {
		 if (e.ff->numVisiblePts() > 0) {
		   auto apex = e.ff->visiblePts(0);
		   if (!apex.isEmpty()) apexes.push_back(apex);
		 }
	       };
    auto fStop = [&]() {return apexes.size() >= n;};
    dfsFacet(H, fVisit, fDo, fStop);
    return apexes;
  }

  facetT* facetWalk() {
    facetT* f = H;
    size_t fSize = f->numVisiblePts();

    auto fVisit = [&](_edge e) {return true;};
    auto fDo = [&](_edge e) {
		 if (e.ff->numVisiblePts() > fSize) {
		   fSize = e.ff->numVisiblePts();
		   f = e.ff;
		 }
	       };
    auto fStop = [&]() { return false; };
    dfsFacet(f, fVisit, fDo, fStop);
    return f;
  }

  /* Choose the furthest outside vertex visible to a facet f
   */
  vertexT furthestApex(facetT *f=nullptr) {
    vertexT apex = vertexT();

    auto fVisit = [&](_edge e) {return true;};
    auto fDo = [&](_edge e) {
		 if (e.ff->numVisiblePts() > 0) apex = e.ff->furthest();
	       };
    auto fStop = [&]() { return !apex.isEmpty(); };
    dfsFacet(f ? f : H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertexT> furthestApexes(size_t n=-1) {
    sequence<vertexT> apexes;

    auto fVisit = [&](_edge e) {return true;};
    auto fDo = [&](_edge e) {
		 if (e.ff->numVisiblePts() > 0) {
		   auto apex = e.ff->furthest();
		   if (!apex.isEmpty()) apexes.push_back(apex);
		 }
	       };
    auto fStop = [&]() {return apexes.size() >= n;};
    dfsFacet(H, fVisit, fDo, fStop);
    return apexes;
  }

  /* Compute a frontier of edges in the clockwise order
      and facets to delete
   */
  tuple<sequence<_edge>*, sequence<facetT*>*> computeFrontier(vertexT apex) {
    facetT* fVisible = apex.attribute.seeFacet;

    auto frontier = new sequence<_edge>();
    auto facets = new sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets->size(); ++i) {
			    if (f == facets->at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](_edge e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || origin.visible(e.fb, apex))
		      return true;
		    else
		      return false;
  		  };

    auto fDo = [&](_edge e) {
		 // Include the facet for deletion if visible
		 bool seeff = origin.visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff))
		   facets->push_back(e.ff);

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 bool seefb = origin.visible(e.fb, apex);
		 if (seefb && !seeff) {
  		   frontier->emplace_back(e.a, e.b, e.ff, e.fb);
		 }
  	       };
    auto fStop = [&](){ return false;};

    dfsEdge(apex.attribute.seeFacet, fVisit, fDo, fStop);
    return make_tuple(frontier, facets);
  }

#ifndef SERIAL
  
  /* Compute a frontier of edges in the clockwise order
      and facets to delete

     Meanwhile reserve the facet using min(apex.attribute.seeFacet*)
   */
  tuple<sequence<_edge>*, sequence<facetT*>*> computeFrontierAndReserve(vertexT apex) {
    facetT* fVisible = apex.attribute.seeFacet;

    auto frontier = new sequence<_edge>();
    auto facets = new sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets->size(); ++i) {
			    if (f == facets->at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](_edge e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || origin.visible(e.fb, apex))
		      return true;
		    else
		      return false;
		  };

    auto fDo = [&](_edge e) {
		 // Include the facet for deletion if visible
		 // Also reserve the facet
		 bool seeff = origin.visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff)) {
		   facets->push_back(e.ff);
		   e.ff->reserve(apex);
		 }

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 // Also reserve invisible facet adjacent to a visible one
		 bool seefb = origin.visible(e.fb, apex);
		 if (seefb && !seeff) {
		   frontier->emplace_back(e.a, e.b, e.ff, e.fb);
		   e.ff->reserve(apex);
		 }
	       };
    auto fStop = [&](){ return false;};

    dfsEdge(apex.attribute.seeFacet, fVisit, fDo, fStop);
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

  bool checkReset() { // todo remove
    bool ok = true;

    auto fVisit = [&](_edge e) { return true; };

    auto fDo = [&](_edge e) {
		 if (e.ff->reservation != -1) ok = false;
	       };

    auto fStop = [&](){ return false; };

    dfsFacet(H, fVisit, fDo, fStop);
    return ok;
  }
#endif // SERIAL

  void redistribute(slice<facetT**, facetT**> facetsBeneath,
		    slice<facetT**, facetT**> newFacets) {

    parlay::write_add(&hSize, newFacets.size() - facetsBeneath.size());

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->numPts();
    }
#ifdef SERIAL
    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->numPts(); ++j) { // Point loop
	facetsBeneath[i]->pts(j).attribute.seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (origin.keep(newFacets[k], facetsBeneath[i]->pts(j))) {
	    facetsBeneath[i]->pts(j).attribute.seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->pts(j), &origin);
	    break;
	  }
	}
      }
    }
#else
    if (fn < 1000) {
      for(int i=0; i<nf; ++i) { // Old facet loop
	for(size_t j=0; j<facetsBeneath[i]->numPts(); ++j) { // Point loop
	  facetsBeneath[i]->pts(j).attribute.seeFacet = nullptr;
	  for (int k=0; k<nnf; ++k) { // New facet loop
	    if (origin.keep(newFacets[k], facetsBeneath[i]->pts(j))) {
	      facetsBeneath[i]->pts(j).attribute.seeFacet = newFacets[k];
	      newFacets[k]->push_back(facetsBeneath[i]->pts(j), &origin);
	      break;
	    }
	  }
	}
      }
    } else {
      auto tmpBuffer = typename facetT::seqT(fn);
      fn = 0;
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
			      if (origin.keep(newFacets[j], tmpBuffer[i])) {
				flag[i] = j;
				tmpBuffer[i].attribute.seeFacet = newFacets[j];
				break;
			      }
			    }
			  });

      auto chunks = split_k(nnf+1, &tmpBuffer, flag);
      for (int j=0; j<nnf; ++j) {
	newFacets[j]->reassign(chunks[j], &origin);
      }
    }
#endif // Serial
  }
  
  std::atomic<size_t>& hullSize() {
    return hSize;
  }

  size_t hullSizeDfs(facetT* start=nullptr) {
    size_t s = 0;
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) { s++;};
    auto fStop = [&]() { return false;};
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    return s;
  }

  // Also checks the hull integrity
  void printHull(facetT* start=nullptr, bool checker=true) {
    if (checker) printHull(start, false);
    size_t hs = hullSizeDfs();
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) {
		 if (checker && hullSizeDfs(e.ff) != hs) {
		   cout << " ..." << endl;
		   cout << "Error, inconsistency detected, abort" << endl;
		   cout << "Erroneous hull size = " << hullSizeDfs(e.ff) << endl;
		   abort();
		 }
		 if (!checker) cout << *e.ff << ":" << e.ff->numVisiblePts() << " ";
	       };
    auto fStop = [&]() { return false;};

    if (!checker) cout << "Hull DFS (" << hs << ")  = ";
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    if (!checker) cout << endl;
  }

  void checkHull(facetT* start=nullptr) {
    cout << "check hull" << endl;
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) {
		 auto f = e.ff;
		 if (!
		     (f->isAdjacent(f->abFacet) && f->abFacet->isAdjacent(f) &&
		      f->isAdjacent(f->bcFacet) && f->bcFacet->isAdjacent(f) &&
		      f->isAdjacent(f->caFacet) && f->caFacet->isAdjacent(f))
		     ) {
		   printHull();
		   throw std::runtime_error("Hull is not linked correctly.");
		 }
	       };
    auto fStop = [&]() { return false;};

    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
  }

  template<class pt>
  void getHull(sequence<facet3d<pt>>& out) {
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) {
		 out.emplace_back(pt((e.ff->a+origin.get()).coords()),
				  pt((e.ff->b+origin.get()).coords()),
				  pt((e.ff->c+origin.get()).coords()));};
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
  }

  template<class pt>
  sequence<pt> getHullPts() {
    sequence<pt> out;
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) {
		 out.emplace_back(pt((e.ff->a+origin.get()).coords()));
		 out.emplace_back(pt((e.ff->b+origin.get()).coords()));
		 out.emplace_back(pt((e.ff->c+origin.get()).coords()));
	       };
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
    parlay::sort_inplace(out);
    return parlay::unique(out);
  }

  // Get all kinds of vertices
  sequence<vertexT> getHullVertices() {
    sequence<vertexT> out;
    auto fVisit = [&](_edge e) { return true;};
    auto fDo = [&](_edge e) {
		 auto f = e.ff;
		 // out.push_back(f->a+origin.get());
		 // out.push_back(f->b+origin.get());
		 // out.push_back(f->c+origin.get());
		 for (auto itr = f->keepList->begin(); itr != f->keepList->end(); itr++) {
		   out.push_back(*itr + origin.get());
		 }
	       };
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
    parlay::sort_inplace(out);
    return parlay::unique(out);
  }

#ifdef WRITE
  void writeHull() {
    using edgeT = _edge;
    ofstream myfile;
    myfile.open("hull.txt", std::ofstream::trunc);
    auto fVisit = [&](edgeT e) { return true;};
    auto fDo = [&](edgeT e) {
		 myfile << e.ff->a << endl;
		 myfile << e.ff->b << endl;
		 myfile << e.ff->c << endl;
	       };
    auto fStop = [&]() { return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    myfile.close();
  }
#endif
};
