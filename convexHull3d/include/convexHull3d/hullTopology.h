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

#include <iostream>
#include <fstream>
#include <atomic>
#include <stack>
#include <tuple>
#include "parlay/sequence.h"
#include "parlay/hash_table.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"

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
class _hullTopology {

protected:

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
      if (!visited(e.ff) && fVisit(e.ff)) {
	fDo(e.ff);
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
      return e1.a==e2.a && e1.b==e2.b;}
  };

  void setHull(facetT* _H) {H = _H;}

/*   facetT* constructorParallel(slice<vertexT*, vertexT*> P) { */

/*     return f0; */
/*   } */

  originT getOrigin() {
    return origin;
  }

 _hullTopology() {}

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

  virtual facetT* initialize(slice<vertexT*, vertexT*> P) = 0;

  std::atomic<size_t>& hullSize() {
    return hSize;
  }

  size_t hullSizeDfs(facetT* start=nullptr) {
    size_t s = 0;
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) { s++;};
    auto fStop = [&]() { return false;};
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    return s;
  }

  // Also checks the hull integrity
  void printHull(facetT* start=nullptr, bool checker=true) {
    if (checker) printHull(start, false);
    size_t hs = hullSizeDfs();
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 if (checker && hullSizeDfs(f) != hs) {
		   cout << " ..." << endl;
		   cout << "Erroneous hull size = " << hullSizeDfs(f) << endl;
		   throw std::runtime_error("Error, hull inconsistency detected");
		 }
		 if (!checker) cout << *f << ":" << f->numVisiblePts() << " ";
	       };
    auto fStop = [&]() { return false;};

    if (!checker) cout << "Hull DFS (" << hs << ")  = ";
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    if (!checker) cout << endl;
  }

  void checkHull(facetT* start=nullptr) {

    auto isAdjacent = [&](facetT* f1, facetT* f2) {
      vertexT V[3]; V[0] = f1->a; V[1] = f2->b; V[2] = f2->c;
      for (int i=0; i<3; ++i) {
	auto v1 = V[i];
	auto v2 = V[(i+1)%3];
	if ( (f2->a == v1 && f2->b == v2) || (f2->a == v2 && f2->b == v1) ) return true;
	if ( (f2->b == v1 && f2->c == v2) || (f2->b == v2 && f2->c == v1) ) return true;
	if ( (f2->c == v1 && f2->a == v2) || (f2->c == v2 && f2->a == v1) ) return true;
      }
      return false;
    };

    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 if (!
		     (isAdjacent(f, f->abFacet) && isAdjacent(f->abFacet, f) &&
		      isAdjacent(f, f->bcFacet) && isAdjacent(f->bcFacet, f) &&
		      isAdjacent(f, f->caFacet) && isAdjacent(f->caFacet, f))
		     ) {
		   printHull();
		   throw std::runtime_error("Hull is not linked correctly.");
		 }
	       };
    auto fStop = [&]() { return false;};

    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
  }

  template<class pt, class facet3d>
  void getHull(sequence<facet3d>& out) {
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 out.emplace_back(pt((f->a+origin.get()).coords()),
				  pt((f->b+origin.get()).coords()),
				  pt((f->c+origin.get()).coords()));};
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
  }

  template<class pt>
  sequence<pt> getHullVertices() {
    sequence<pt> out;
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 out.emplace_back(pt((f->a+origin.get()).coords()));
		 out.emplace_back(pt((f->b+origin.get()).coords()));
		 out.emplace_back(pt((f->c+origin.get()).coords()));
  	       };
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
    parlay::sort_inplace(out);
    return parlay::unique(out);
  }

  sequence<vertexT> getHullPts() {
    sequence<vertexT> out;
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 for (size_t i = 0; i < f->numPts(); ++i) {
		   out.push_back(f->pts(i) + origin.get());
		 }
	       };
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
    parlay::sort_inplace(out);
    return parlay::unique(out);
  }

  void writeHull(char const *fileName) {
    using edgeT = _edge;
    auto offset = origin.get();
    ofstream myfile;
    myfile.open(fileName, std::ofstream::trunc);
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 myfile << f->a + offset << endl;
		 myfile << f->b + offset << endl;
		 myfile << f->c + offset << endl;
	       };
    auto fStop = [&]() { return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    myfile.close();
  }
};
