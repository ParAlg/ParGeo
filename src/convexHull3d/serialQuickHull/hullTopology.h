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

#include <atomic>
#include <fstream>
#include <iostream>
#include <stack>
#include <tuple>
#include <vector>
#include <unordered_map>
#include "pargeo/getTime.h"
#include "pargeo/point.h"

namespace pargeo::hull3d {

template <class facetT, class vertexT>
class _hullTopology {

protected:

  static constexpr typename vertexT::floatT eps = vertexT::pointT::eps;

public:

  // The number of facets in H
  std::atomic<size_t> hSize;

  vertexT interiorPt;

  // A linked structure for the facets
  facetT* H;

  inline bool visible(facetT* f, vertexT p) {
    return (f->a - p).dot(f->area) > eps;
  }

  /* Depth-first hull traversal (no facet repeat)
  */
  template <class F, class G, class H>
  void dfsFacet(facetT* start, F& fVisit, G& fDo, H& fStop) {

    std::unordered_map<facetT*, int> V;
    auto mark = [&](facetT* f) { V[f] = 1; };
    auto visited = [&](facetT* f) { return V.find(f) != V.end(); };

    std::stack<_edge> S;

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
    std::vector<edgeT> V;
    auto mark = [&](edgeT f) {V.push_back(f);};
    auto visited = [&](edgeT f) {
      for (auto g: V) {
	if (f == g) return true;
      }
      return false;};

    std::stack<edgeT> S;

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

  // An arbitrary coordinate class located within the hull
  // it also contains primitives for visibility test
  // originT origin;

  // todo test this function
  // try to link facets f1 and f2, and return success
  bool linkFacet(facetT* f1, facetT* f2) {
    bool matched = false;

    auto match = [&](vertexT f1v1,
		     vertexT f1v2,
		     vertexT f2v1,
		     vertexT f2v2) {
		   if ((f1v1 == f2v1 && f1v2 == f2v2) ||
		       (f1v1 == f2v2 && f1v2 == f2v1)) {
		     matched = true;
		     return true;
		   }
		   return false;
		 };

    if (!matched && match(f1->a, f1->b, f2->a, f2->b)) {
      f1->abFacet = f2; f2->abFacet = f1;}
    if (!matched && match(f1->b, f1->c, f2->a, f2->b)) {
      f1->bcFacet = f2; f2->abFacet = f1;}
    if (!matched && match(f1->c, f1->a, f2->a, f2->b)) {
      f1->caFacet = f2; f2->abFacet = f1;}

    if (!matched && match(f1->a, f1->b, f2->b, f2->c)) {
      f1->abFacet = f2; f2->bcFacet = f1;}
    if (!matched && match(f1->b, f1->c, f2->b, f2->c)) {
      f1->bcFacet = f2; f2->bcFacet = f1;}
    if (!matched && match(f1->c, f1->a, f2->b, f2->c)) {
      f1->caFacet = f2; f2->bcFacet = f1;}

    if (!matched && match(f1->a, f1->b, f2->c, f2->a)) {
      f1->abFacet = f2; f2->caFacet = f1;}
    if (!matched && match(f1->b, f1->c, f2->c, f2->a)) {
      f1->bcFacet = f2; f2->bcFacet = f1;}
    if (!matched && match(f1->c, f1->a, f2->c, f2->a)) {
      f1->caFacet = f2; f2->caFacet = f1;}

    return matched;
  }

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

  _hullTopology() {}

  /* Compute a frontier of edges in the clockwise order
      and facets to delete
   */
  std::tuple<std::vector<_edge>, std::vector<facetT*>> computeFrontier(vertexT apex) {

    facetT* fVisible = apex.seeFacet;

    auto frontier = std::vector<_edge>();
    auto facets = std::vector<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets.size(); ++i) {
			    if (f == facets.at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](_edge e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || visible(e.fb, apex))
		      return true;
		    else
		      return false;
  		  };

    auto fDo = [&](_edge e) {
		 // Include the facet for deletion if visible
		 bool seeff = visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff))
		   facets.push_back(e.ff);

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 bool seefb = visible(e.fb, apex);
		 if (seefb && !seeff) {
		   frontier.emplace_back(e.a, e.b, e.ff, e.fb);
		 }
  	       };
    auto fStop = [&](){ return false;};

    dfsEdge(apex.seeFacet, fVisit, fDo, fStop);
    return std::make_tuple(std::move(frontier), std::move(facets));
  }

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
		   std::cout << " ...\n";
		   std::cout << "Erroneous hull size = " << hullSizeDfs(f) << "\n";
		   throw std::runtime_error("Error, hull inconsistency detected");
		 }
		 if (!checker) std::cout << f << ":" << f->numPts() << " ";
	       };
    auto fStop = [&]() { return false;};

    if (!checker) std::cout << "Hull DFS (" << hs << ")  = ";
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    if (!checker) std::cout << "\n";
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

  template <class seqT>
  void getFacet(seqT& out) {
    using pt = vertexT;
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 out.emplace_back(pt((f->a+interiorPt).coords()),
				  pt((f->b+interiorPt).coords()),
				  pt((f->c+interiorPt).coords()));};
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
  }

  void writeHull(char const *fileName) {
    using edgeT = _edge;
    std::ofstream myfile;
    myfile.open(fileName, std::ofstream::trunc);
    auto fVisit = [&](facetT* f) { return true;};
    auto fDo = [&](facetT* f) {
		 myfile << f->a + interiorPt << "\n";
		 myfile << f->b + interiorPt << "\n";
		 myfile << f->c + interiorPt << "\n";
	       };
    auto fStop = [&]() { return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    myfile.close();
  }
};

} // End namespaces
