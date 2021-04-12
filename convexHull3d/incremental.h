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

//#define WRITE // Write to file, visualize using python3 plot.py
//#define VERBOSE
#define SILENT

#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

#include <set>
#include <atomic>
#include <vector>
#include <assert.h>
#include <stack>
#include <tuple>
#include "parlay/sequence.h"
#include "parlay/hash_table.h"
#include "geometry/point.h"
#include "geometry/algebra.h"
#include "common/get_time.h"
#include "hull.h"
#include "split.h"
#include "shared.h"

using namespace parlay;
using namespace parlay::internal;
using pointT = pargeo::fpoint<3>;
static const typename pointT::floatT numericKnob = 1e-5;

////////////////////////////
// Facet and vertex with attributes
////////////////////////////

struct linkedFacet3d;

struct vertexAtt;

using vertex3d = pargeo::_point<3, pointT::floatT, pointT::floatT, vertexAtt>;

struct vertexAtt {
#ifdef VERBOSE
  size_t i;
#endif
  linkedFacet3d *seeFacet;
  vertexAtt() {}
};

static std::ostream& operator<<(std::ostream& os, const vertex3d& v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

/* Signed volume (x6) of an oriented tetrahedron (example below is positive).
     d
     o
    /|\
   / | o b
  o--o/
   a   c

    x
     origin
*/

template <class pt>
inline typename pt::floatT signedVolume(pt a, pt b, pt c, pt d) {
  return (a-d).dot(crossProduct3d(b-a, c-a));
}

// area is precomputed from oriented triang a,b,c
template <class pt>
inline typename pt::floatT signedVolume(pt a, pt d, pt area) {
  return (a-d).dot(area);
}

template<class facetT, class vertexT>
bool visible(facetT* f, vertexT p) {
  // if (signedVolume(f->a, f->b, f->c, p) > numericKnob)
  if (signedVolume(f->a, p, f->area) > numericKnob)
    return true;
  else
    return false;
}

template<class facetT, class vertexT>
bool visibleNoDup(facetT* f, vertexT p) {
  if (signedVolume(f->a, p, f->area) > numericKnob)
    return true && f->a != p && f->b != p && f->c != p;
  else
    return false;
}

// todo check if still necessary
template<class facetT, class vertexT>
bool visibleCast(facetT* f, vertexT p) {
  if (signedVolume(vertex3d(f->a.coords()),
		   vertex3d(f->b.coords()),
		   vertex3d(f->c.coords()),
		   vertex3d(p.coords())
		   ) > numericKnob)
    return true;
  else
    return false;
}

struct linkedFacet3d {
#ifdef SERIAL
  typedef vector<vertex3d> seqT;
#else
  typedef sequence<vertex3d> seqT;
#endif

  vertex3d a, b, c;
  linkedFacet3d *abFacet;
  linkedFacet3d *bcFacet;
  linkedFacet3d *caFacet;
  vertex3d area;

  bool isAdjacent(linkedFacet3d* f2) {
    vertex3d V[3]; V[0] = a; V[1] = b; V[2] = c;
    for (int i=0; i<3; ++i) {
      auto v1 = V[i];
      auto v2 = V[(i+1)%3];
      if ( (f2->a == v1 && f2->b == v2) || (f2->a == v2 && f2->b == v1) ) return true;
      if ( (f2->b == v1 && f2->c == v2) || (f2->b == v2 && f2->c == v1) ) return true;
      if ( (f2->c == v1 && f2->a == v2) || (f2->c == v2 && f2->a == v1) ) return true;
    }
    return false;
  }

  // Stores the minimum memory address of the seeFacet of the reserving vertices
  std::atomic<size_t> reservation; // todo (not always used but has little impact on speed)

  void reserve(vertex3d& p) {
    parlay::write_min(&reservation, (size_t)p.attribute.seeFacet, std::less<size_t>());
  }

  bool reserved(vertex3d& p) {
    return reservation == (size_t)p.attribute.seeFacet;
  }

  seqT *seeList;

  vertex3d& at(size_t i) { return seeList->at(i); }

  size_t size() { return seeList->size(); }

  void clear() { seeList->clear(); }

  void push_back(vertex3d v) { seeList->push_back(v); }

  vertex3d furthest() {
    typename vertex3d::floatT m = numericKnob;
    auto apex = vertex3d();

#ifdef SERIAL
    for (size_t i=0; i<size(); ++i) {
      auto m2 = signedVolume(a, b, c, at(i));
      if (m2 > m) {
	m = m2;
	apex = at(i);
      }
    }
#else
    apex = parlay::max_element(seeList->cut(0, seeList->size()),
			       [&](vertex3d aa, vertex3d bb) {
				 return signedVolume(a, b, c, aa) < signedVolume(a, b, c, bb);
			       });
#endif
    return apex;
  }

  linkedFacet3d(vertex3d _a, vertex3d _b, vertex3d _c): a(_a), b(_b), c(_c) {
    if (pargeo::determinant3by3(a, b, c) > numericKnob)
      swap(b, c);

    seeList = new seqT();

    reservation = -1; // (unsigned)

    area = crossProduct3d(b-a, c-a);
  }

  ~linkedFacet3d() {
    delete seeList;
  }
};

static std::ostream& operator<<(std::ostream& os, const linkedFacet3d& v) {
#ifdef VERBOSE
  os << "(" << v.a.attribute.i << "," << v.b.attribute.i << "," << v.c.attribute.i << ")";
#else
  os << "(" << v.a << "," << v.b << "," << v.c << ")";
#endif
  return os;
}

// a, b, c can be input in any order
linkedFacet3d* makeFacet(vertex3d a, vertex3d b, vertex3d c) {
  auto f = new linkedFacet3d(a, b, c);
  return f;
}

// Link f with ab, bc, ca; the edge matching is automatic -- input in any order
template<class fc>
void linkFacet(fc* f, fc* ab, fc* bc, fc* ca) {
  fc* F[3]; F[0]=ab; F[1]=bc; F[2]=ca;
  auto findFacet = [&](vertex3d v1, vertex3d v2) {
		     for(int i=0; i<3; ++i) {
		       if ((F[i]->a==v1 && F[i]->b==v2) || (F[i]->b==v1 && F[i]->a==v2) ||
			   (F[i]->b==v1 && F[i]->c==v2) || (F[i]->c==v1 && F[i]->b==v2) ||
			   (F[i]->c==v1 && F[i]->a==v2) || (F[i]->a==v1 && F[i]->c==v2)) {
			 //cout << "linking " << *F[i] << " to f's (" << v1.attribute.i << "," << v2.attribute.i << ")" << endl;
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

  auto linkEdge = [&](fc* f1, fc* f2, vertex3d v1, vertex3d v2) {
		    if ( (f2->a==v1 && f2->b==v2) || (f2->a==v2 && f2->b==v1) ) {
		      // cout << "linking back (" << f2->a.attribute.i << " " << f2->b.attribute.i << ") to " << *f1 << endl;
		      f2->abFacet = f1;
		    } else if ( (f2->b==v1 && f2->c==v2) || (f2->b==v2 && f2->c==v1) ) {
		      // cout << "linking back (" << f2->b.attribute.i << " " << f2->c.attribute.i << ") to " << *f1 << endl;
		      f2->bcFacet = f1;
		    } else if ( (f2->c==v1 && f2->a==v2) || (f2->c==v2 && f2->a==v1) ) {
		      // cout << "linking back (" << f2->c.attribute.i << " " << f2->a.attribute.i << ") to " << *f1 << endl;
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

////////////////////////////
// Hull structure
////////////////////////////

/* An edge structure for advancing direction in a traversal

   ff: forward facet to advance to
   fb: parent facet from which we come

       o
      / \ --->ff
   a o---o b
      \ / --->fb
       o
*/
template <class facetT, class vertexT>
struct _edge {
  vertexT a, b;
  facetT* ff;
  facetT* fb;

  _edge(vertexT _a, vertexT _b, facetT* _ff, facetT* _fb) {
    a = _a; b = _b; ff = _ff; fb = _fb;}

  friend bool operator==(_edge e1, _edge e2) {
    return e1.a==e2.a && e1.b==e2.b;}//todo
};

/* A structure that maintains and traverses the facets
 */
template <class facetT, class vertexT>
struct _hull {
  vertexT origin;
  facetT* H;
  std::atomic<size_t> hSize;

  // Initialize from a given linked-facet-hull and an origin
  template<class pt>
  _hull(facetT* _H, pt _origin, size_t _hSize) {
    H = _H;
    hSize = _hSize;
    for (int i=0; i<3; ++i)
      origin[i] = _origin[i];
  }

  // Initialize from a point array and four vertices of a simplex
  template<class pt>
  _hull(slice<pt*, pt*> P, size_t c1, size_t c2, size_t c3, size_t c4) {
    hSize = 4;

    pt tmp = P[c1] + P[c2] + P[c3] + P[c4];
    origin[0] = tmp[0] / 4;
    origin[1] = tmp[1] / 4;
    origin[2] = tmp[2] / 4;

    // Initialize points with visible facet link
    auto Q = linkedFacet3d::seqT(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				for(int j=0; j<P[i].dim; ++j)
				  Q[i][j] = P[i][j] - origin[j];
#ifdef VERBOSE
				Q[i].attribute.i = i;
#endif
			      });

    // Make initial facets
    auto f0 = makeFacet(Q[c1], Q[c2], Q[c3]);
    auto f1 = makeFacet(Q[c1], Q[c2], Q[c4]);
    auto f2 = makeFacet(Q[c3], Q[c4], Q[c2]);
    auto f3 = makeFacet(Q[c3], Q[c4], Q[c1]);

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
	if (visibleNoDup(f0, Q[i])) {
	  Q[i].attribute.seeFacet = f0;
	  f0->push_back(Q[i]);
	} else if (visibleNoDup(f1, Q[i])) {
	  Q[i].attribute.seeFacet = f1;
	  f1->push_back(Q[i]);
	} else if (visibleNoDup(f2, Q[i])) {
	  Q[i].attribute.seeFacet = f2;
	  f2->push_back(Q[i]);
	} else if (visibleNoDup(f3, Q[i])) {
	  Q[i].attribute.seeFacet = f3;
	  f3->push_back(Q[i]);
	} else {
	  Q[i].attribute.seeFacet = nullptr;
	}
      }

    } else {

      auto flag = sequence<int>(P.size());
      parallel_for(0, P.size(), [&](size_t i) {
				  if (visibleNoDup(f0, Q[i])) {
				    flag[i] = 0; Q[i].attribute.seeFacet = f0;
				  } else if (visibleNoDup(f1, Q[i])) {
				    flag[i] = 1; Q[i].attribute.seeFacet = f1;
				  } else if (visibleNoDup(f2, Q[i])) {
				    flag[i] = 2; Q[i].attribute.seeFacet = f2;
				  } else if (visibleNoDup(f3, Q[i])) {
				    flag[i] = 3; Q[i].attribute.seeFacet = f3;
				  } else {
				    flag[i] = 4; Q[i].attribute.seeFacet = nullptr;
				  }
				});

      auto chunks = split_k(5, &Q, flag);

      f0->seeList = chunks[0];
      f1->seeList = chunks[1];
      f2->seeList = chunks[2];
      f3->seeList = chunks[3];
    }

#ifdef VERBOSE
    cout << "initial-hull:" << endl;
    cout << " facet 0, area = " << f0->area.length()/2 << " #pts = " << f0->size() << endl;
    cout << " facet 1, area = " << f1->area.length()/2 << " #pts = " << f1->size() << endl;
    cout << " facet 2, area = " << f2->area.length()/2 << " #pts = " << f2->size() << endl;
    cout << " facet 3, area = " << f3->area.length()/2 << " #pts = " << f3->size() << endl;
    cout << " volume = " << abs((f0->a-Q[c4]).dot(f0->area)/6) << endl;
    // cout << "vertices: " << endl;
    // cout << " v0 = " << Q[c1]+origin << endl;
    // cout << " v1 = " << Q[c2]+origin << endl;
    // cout << " v2 = " << Q[c3]+origin << endl;
    // cout << " v3 = " << Q[c4]+origin << endl;
    // cout << "origin = " << origin << endl;
#endif

    H = f0;
  }

  /* Depth-first hull traversal (no facet repeat)
  */
  template <class F, class G, class H>
  void dfsFacet(facetT* start, F& fVisit, G& fDo, H& fStop) {
    // sequence<facetT*> V;
    // auto mark = [&](facetT* f) {V.push_back(f);};
    // auto visited = [&](facetT* f) {
    // 		     for (auto g: V) {
    // 		       if (f == g) return true;
    // 		     }
    // 		     return false;};

    hashtable<hash_pointer<facetT*>> V(hSize, hash_pointer<facetT*>());
    auto mark = [&](facetT* f) {V.insert(f);};
    auto visited = [&](facetT* f) {return V.find(f) != nullptr;};

    stack<_edge<facetT, vertexT>> S;

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
    using edgeT = _edge<facetT, vertexT>;

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

  std::atomic<size_t>& hullSize() {
    return hSize;
  }

  size_t hullSizeDfs(facetT* start=nullptr) {
    size_t s = 0;
    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) { s++;};
    auto fStop = [&]() { return false;};
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    return s;
  }

  // Also checks the hull integrity
  void printHull(facetT* start=nullptr, bool checker=true) {
    if (checker) printHull(start, false);
    size_t hs = hullSizeDfs();
    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (checker && hullSizeDfs(e.ff) != hs) {
		   cout << " ..." << endl;
		   cout << "Error, inconsistency detected, abort" << endl;
		   cout << "Erroneous hull size = " << hullSizeDfs(e.ff) << endl;
		   abort();
		 }
		 if (!checker) cout << *e.ff << ":" << e.ff->size() << " ";
	       };
    auto fStop = [&]() { return false;};

    if (!checker) cout << "Hull DFS (" << hs << ")  = ";
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    if (!checker) cout << endl;
  }

  void checkHull(facetT* start=nullptr) {
    cout << "check hull" << endl;
    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
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

    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 out.emplace_back(pt((e.ff->a+origin).coords()),
				  pt((e.ff->b+origin).coords()),
				  pt((e.ff->c+origin).coords()));};
    auto fStop = [&]() { return false;};

    dfsFacet(H, fVisit, fDo, fStop);
  }

#ifdef WRITE
  void writeHull() {
    using edgeT = _edge<facetT, vertexT>;
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

////////////////////////////
// Conflict graph
////////////////////////////

template <class facetT, class vertexT>
class conflictList {

public:
  _hull<facetT, vertexT> *context;

  conflictList(_hull<facetT, vertexT> *_c): context(_c) {};

  /* Choose a random outside vertex visible to some facet
   */
  vertexT randomApex() {
    if (context->H->size() > 0) return context->H->at(0);
    vertexT apex = vertexT();
    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->size() > 0)
		   apex = e.ff->at(0);};
    auto fStop = [&]() {
		   if (!apex.isEmpty()) return true;
		   else return false;};
    context->dfsFacet(context->H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertex3d> randomApexes(size_t n=-1) {
    sequence<vertex3d> apexes;

    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->size() > 0) {
		   auto apex = e.ff->at(0);
		   if (!apex.isEmpty()) apexes.push_back(apex);
		 }
	       };
    auto fStop = [&]() {return apexes.size() >= n;};
    context->dfsFacet(context->H, fVisit, fDo, fStop);
    return apexes;
  }

  facetT* facetWalk() {
    facetT* f = context->H;
    size_t fSize = f->size();

    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->size() > fSize) {
		   fSize = e.ff->size();
		   f = e.ff;
		 }
	       };
    auto fStop = [&]() { return false; };
    context->dfsFacet(f, fVisit, fDo, fStop);
    return f;
  }

  /* Choose the furthest outside vertex visible to a facet f
   */
  vertexT furthestApex(facetT *f=nullptr) {
    vertexT apex = vertexT();

    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->size() > 0) apex = e.ff->furthest();
	       };
    auto fStop = [&]() { return !apex.isEmpty(); };
    context->dfsFacet(f ? f : context->H, fVisit, fDo, fStop);
    return apex;
  }

  // If n not supplied, find furthest apexes of all facets
  sequence<vertex3d> furthestApexes(size_t n=-1) {
    sequence<vertex3d> apexes;

    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->size() > 0) {
		   auto apex = e.ff->furthest();
		   if (!apex.isEmpty()) apexes.push_back(apex);
		 }
	       };
    auto fStop = [&]() {return apexes.size() >= n;};
    context->dfsFacet(context->H, fVisit, fDo, fStop);
    return apexes;
  }

  /* Compute a frontier of edges in the clockwise order
      and facets to delete
   */
  tuple<sequence<_edge<facetT, vertexT>>*, sequence<facetT*>*> computeFrontier(vertexT apex) {
    facetT* fVisible = apex.attribute.seeFacet;

    auto frontier = new sequence<_edge<facetT, vertexT>>();
    auto facets = new sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets->size(); ++i) {
			    if (f == facets->at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](_edge<facetT, vertexT> e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || visible(e.fb, apex))
		      return true;
		    else
		      return false;
  		  };

    auto fDo = [&](_edge<facetT, vertexT> e) {
		 // Include the facet for deletion if visible
		 bool seeff = visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff))
		   facets->push_back(e.ff);

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 bool seefb = visible(e.fb, apex);
		 if (seefb && !seeff) {
  		   frontier->emplace_back(e.a, e.b, e.ff, e.fb);
		 }
  	       };
    auto fStop = [&](){ return false;};

    context->dfsEdge(apex.attribute.seeFacet, fVisit, fDo, fStop);
    return make_tuple(frontier, facets);
  }

  /* Compute a frontier of edges in the clockwise order
      and facets to delete

     Meanwhile reserve the facet using min(apex.attribute.seeFacet*)
   */
  tuple<sequence<_edge<facetT, vertexT>>*, sequence<facetT*>*> computeFrontierAndReserve(vertexT apex) {
    facetT* fVisible = apex.attribute.seeFacet;

    auto frontier = new sequence<_edge<facetT, vertexT>>();
    auto facets = new sequence<facetT*>();
    auto facetVisited = [&](facetT* f) {
			  for (size_t i=0; i<facets->size(); ++i) {
			    if (f == facets->at(i)) return true;
			  }
			  return false;
			};

    auto fVisit = [&](_edge<facetT, vertexT> e) {
		    // Visit the facet as long as the parent facet is visible to the apex
		    // e.fb == nullptr for the starting facet (whose parent is nullptr, see dfsEdge(...))
		    if (e.fb == nullptr || visible(e.fb, apex))
		      return true;
		    else
		      return false;
		  };

    auto fDo = [&](_edge<facetT, vertexT> e) {
		 // Include the facet for deletion if visible
		 // Also reserve the facet
		 bool seeff = visible(e.ff, apex);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff)) {
		   facets->push_back(e.ff);
		   e.ff->reserve(apex);
		 }

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 // Also reserve invisible facet adjacent to a visible one
		 bool seefb = visible(e.fb, apex);
		 if (seefb && !seeff) {
		   frontier->emplace_back(e.a, e.b, e.ff, e.fb);
		   e.ff->reserve(apex);
		 }
	       };
    auto fStop = [&](){ return false;};

    context->dfsEdge(apex.attribute.seeFacet, fVisit, fDo, fStop);
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

    auto fVisit = [&](_edge<facetT, vertexT> e) { return true; };

    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->reservation != -1) ok = false;
	       };

    auto fStop = [&](){ return false; };

    context->dfsFacet(context->H, fVisit, fDo, fStop);
    return ok;
  }

  void redistribute(slice<facetT**, facetT**> facetsBeneath,
		    slice<facetT**, facetT**> newFacets) {

    parlay::write_add(&context->hSize, newFacets.size() - facetsBeneath.size());

    // Redistribute the outside points

    int nf = facetsBeneath.size();
    int nnf = newFacets.size();

    size_t fn = 0;
    for(int j=0; j<nf; ++j) {
      fn += facetsBeneath[j]->size();
    }
#ifdef SERIAL
    for(int i=0; i<nf; ++i) { // Old facet loop
      for(size_t j=0; j<facetsBeneath[i]->size(); ++j) { // Point loop
	facetsBeneath[i]->at(j).attribute.seeFacet = nullptr;
	for (int k=0; k<nnf; ++k) { // New facet loop
	  if (visibleNoDup(newFacets[k], facetsBeneath[i]->at(j))) {
	    facetsBeneath[i]->at(j).attribute.seeFacet = newFacets[k];
	    newFacets[k]->push_back(facetsBeneath[i]->at(j));
	    break;
	  }
	}
      }
    }
#else
    if (fn < 1000) {
      for(int i=0; i<nf; ++i) { // Old facet loop
	for(size_t j=0; j<facetsBeneath[i]->size(); ++j) { // Point loop
	  facetsBeneath[i]->at(j).attribute.seeFacet = nullptr;
	  for (int k=0; k<nnf; ++k) { // New facet loop
	    if (visibleNoDup(newFacets[k], facetsBeneath[i]->at(j))) {
	      facetsBeneath[i]->at(j).attribute.seeFacet = newFacets[k];
	      newFacets[k]->push_back(facetsBeneath[i]->at(j));
	      break;
	    }
	  }
	}
      }
    } else {
      auto tmpBuffer = typename linkedFacet3d::seqT(fn);
      fn = 0;
      for(int j=0; j<nf; ++j) {
	parallel_for(0, facetsBeneath[j]->size(),
		     [&](size_t x){tmpBuffer[fn+x] = facetsBeneath[j]->at(x);});
	fn += facetsBeneath[j]->size();
      }

      auto flag = sequence<int>(fn);
      parallel_for(0, fn, [&](size_t i) {
			    flag[i] = nnf;
			    tmpBuffer[i].attribute.seeFacet = nullptr;
			    for (int j=0; j<nnf; ++j) {
			      if (visibleNoDup(newFacets[j], tmpBuffer[i])) {
				flag[i] = j;
				tmpBuffer[i].attribute.seeFacet = newFacets[j];
				break;
			      }
			    }
			  });

      auto chunks = split_k(nnf+1, &tmpBuffer, flag);
      for (int j=0; j<nnf; ++j) {
	delete newFacets[j]->seeList;
	newFacets[j]->seeList = chunks[j];
      }
    }
#endif // Serial

  }
};

////////////////////////////
// Initialization
////////////////////////////

template<class pt>
conflictList<linkedFacet3d, vertex3d> *makeInitial8Hull(slice<pt*, pt*> P) {

  // Find 8 extrema assuming no duplicate
  auto xx = minmax_element(P, [&](pt i, pt j) {return i[0]<j[0];});
  size_t xMin = xx.first - &P[0];
  size_t xMax = xx.second - &P[0];
  assert(xMin != xMax);

  auto yy = minmax_element(P, [&](pt i, pt j) {return i[1]<j[1];});
  size_t yMin = yy.first - &P[0];
  size_t yMax = yy.second - &P[0];
  assert(yMin != yMax);

  auto zz = minmax_element(P, [&](pt i, pt j) {return i[2]<j[2];});
  size_t zMin = zz.first - &P[0];
  size_t zMax = zz.second - &P[0];
  assert(zMin != zMax);

  pt origin = P[xMin] + P[xMax] + P[yMin] + P[yMax] + P[zMin] + P[zMax];
  origin[0] = origin[0] / 6;
  origin[1] = origin[1] / 6;
  origin[2] = origin[2] / 6;

  // Initialize points with visible facet link
  auto Q = linkedFacet3d::seqT(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      for(int j=0; j<P[i].dim; ++j)
				Q[i][j] = P[i][j] - origin[j];
#ifdef VERBOSE
			      Q[i].attribute.i = i;
#endif
			    });

  // Make initial facets
  auto f0 = makeFacet(Q[xMin], Q[yMin], Q[zMin]);
  auto f1 = makeFacet(Q[xMax], Q[yMin], Q[zMin]);
  auto f2 = makeFacet(Q[xMin], Q[yMax], Q[zMin]);
  auto f3 = makeFacet(Q[xMax], Q[yMax], Q[zMin]);
  auto f4 = makeFacet(Q[xMin], Q[yMin], Q[zMax]);
  auto f5 = makeFacet(Q[xMax], Q[yMin], Q[zMax]);
  auto f6 = makeFacet(Q[xMin], Q[yMax], Q[zMax]);
  auto f7 = makeFacet(Q[xMax], Q[yMax], Q[zMax]);

  linkFacet(f0, f2, f1, f4);
  linkFacet(f1, f5, f0, f3);
  linkFacet(f2, f6, f3, f0);
  linkFacet(f3, f1, f2, f7);
  linkFacet(f4, f0, f5, f6);
  linkFacet(f5, f7, f4, f1);
  linkFacet(f6, f4, f7, f2);
  linkFacet(f7, f3, f6, f5);

#ifdef SERIAL
  bool serial = true;
#else
  bool serial = false;
#endif

  if (Q.size() < 1000 || serial) {

    for(size_t i=0; i<P.size(); i++) {
      if (visibleNoDup(f0, Q[i])) {
	Q[i].attribute.seeFacet = f0;
	f0->push_back(Q[i]);
      } else if (visibleNoDup(f1, Q[i])) {
	Q[i].attribute.seeFacet = f1;
	f1->push_back(Q[i]);
      } else if (visibleNoDup(f2, Q[i])) {
	Q[i].attribute.seeFacet = f2;
	f2->push_back(Q[i]);
      } else if (visibleNoDup(f3, Q[i])) {
	Q[i].attribute.seeFacet = f3;
	f3->push_back(Q[i]);
      } else if (visibleNoDup(f4, Q[i])) {
	Q[i].attribute.seeFacet = f4;
	f4->push_back(Q[i]);
      } else if (visibleNoDup(f5, Q[i])) {
	Q[i].attribute.seeFacet = f5;
	f5->push_back(Q[i]);
      } else if (visibleNoDup(f6, Q[i])) {
	Q[i].attribute.seeFacet = f6;
	f6->push_back(Q[i]);
      } else if (visibleNoDup(f7, Q[i])) {
	Q[i].attribute.seeFacet = f7;
	f7->push_back(Q[i]);
      } else {
	Q[i].attribute.seeFacet = nullptr;
      }
    }

  } else {

    auto flag = sequence<int>(P.size());
    parallel_for(0, P.size(), [&](size_t i) {
				if (visibleNoDup(f0, Q[i])) {
				  flag[i] = 0; Q[i].attribute.seeFacet = f0;
				} else if (visibleNoDup(f1, Q[i])) {
				  flag[i] = 1; Q[i].attribute.seeFacet = f1;
				} else if (visibleNoDup(f2, Q[i])) {
				  flag[i] = 2; Q[i].attribute.seeFacet = f2;
				} else if (visibleNoDup(f3, Q[i])) {
				  flag[i] = 3; Q[i].attribute.seeFacet = f3;
				} else if (visibleNoDup(f4, Q[i])) {
				  flag[i] = 4; Q[i].attribute.seeFacet = f4;
				} else if (visibleNoDup(f5, Q[i])) {
				  flag[i] = 5; Q[i].attribute.seeFacet = f5;
				} else if (visibleNoDup(f6, Q[i])) {
				  flag[i] = 6; Q[i].attribute.seeFacet = f6;
				} else if (visibleNoDup(f7, Q[i])) {
				  flag[i] = 7; Q[i].attribute.seeFacet = f7;
				} else {
				  flag[i] = 8; Q[i].attribute.seeFacet = nullptr;
				}
			      });

    auto chunks = split_k(9, &Q, flag);

    f0->seeList = chunks[0];
    f1->seeList = chunks[1];
    f2->seeList = chunks[2];
    f3->seeList = chunks[3];
    f4->seeList = chunks[4];
    f5->seeList = chunks[5];
    f6->seeList = chunks[6];
    f7->seeList = chunks[7];
  }

  auto context = new _hull<linkedFacet3d, vertex3d>(f0, origin);

  return new conflictList(context);
}

template<class pt>
conflictList<linkedFacet3d, vertex3d> *makeInitialSimplex0(slice<pt*, pt*> P) {

  size_t X[6]; // extrema

  auto xx = minmax_element(P, [&](pt i, pt j) {return i[0]<j[0];});
  X[0] = xx.first - &P[0]; X[1] = xx.second - &P[0];

  auto yy = minmax_element(P, [&](pt i, pt j) {return i[1]<j[1];});
  X[2] = yy.first - &P[0]; X[3] = yy.second - &P[0];

  auto zz = minmax_element(P, [&](pt i, pt j) {return i[2]<j[2];});
  X[4] = zz.first - &P[0]; X[5] = zz.second - &P[0];

  set<size_t, greater<size_t>> pts;
  if (P[X[1]][0]-P[X[0]][0] > P[X[3]][1]-P[X[2]][1] && P[X[1]][0]-P[X[0]][0] > P[X[5]][2]-P[X[4]][2]) {
    pts.insert(X[0]); pts.insert(X[1]); pts.insert(X[3]); pts.insert(X[5]);
  } else if (P[X[3]][1]-P[X[2]][1] > P[X[1]][0]-P[X[0]][0] && P[X[3]][1]-P[X[2]][1] > P[X[5]][2]-P[X[4]][2]) {
    pts.insert(X[2]); pts.insert(X[3]); pts.insert(X[1]); pts.insert(X[5]);
  } else {
    pts.insert(X[4]); pts.insert(X[5]); pts.insert(X[1]); pts.insert(X[3]);
  }

  int pp = 0;
  while (pts.size() < 4) pts.insert(X[pp++]);

  auto itr = pts.begin();
  size_t c1 = *itr; itr++; size_t c2 = *itr; itr++;
  size_t c3 = *itr; itr++; size_t c4 = *itr;

  auto context = new _hull<linkedFacet3d, vertex3d>(P, c1, c2, c3, c4);

  return new conflictList(context);
}

template<class pt>
conflictList<linkedFacet3d, vertex3d> *makeInitialSimplex1(slice<pt*, pt*> P) {

  // Pick largest triangle from sample
  srand(1);
  typename pt::floatT myArea = 0;
  auto keep = [&](pt a, pt b, pt c) {
		typename pt::floatT newArea = crossProduct3d(b-a, c-a).length();
		if (newArea > myArea) {
		  myArea = newArea;
		  return true;
		} else
		  return false;
	      };
  size_t xMin, xMax, yApex;
  pt x1, x2, y1;
  for (size_t i=0; i<P.size()/100; ++i) {
    size_t xMin2 = rand()%P.size();
    size_t xMax2 = rand()%P.size();
    size_t yApex2 = rand()%P.size();
    if (keep(P[xMin2], P[xMax2], P[yApex2])) {
      xMin = xMin2;
      xMax = xMax2;
      yApex = yApex2;
    }
  }
  x1 = P[xMin];
  x2 = P[xMax];
  y1 = P[yApex];

  // Maximize simplex volume
  pt area = crossProduct3d(x1-y1, x2-y1);
  auto z = max_element(P, [&](pt i, pt j) {
			    return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
			  });
  size_t zApex = z - &P[0];

  size_t c1 = xMin;
  size_t c2 = xMax;
  size_t c3 = yApex;
  size_t c4 = zApex;

  auto context = new _hull<linkedFacet3d, vertex3d>(P, c1, c2, c3, c4);

  return new conflictList(context);
}

template<class pt>
conflictList<linkedFacet3d, vertex3d> *makeInitialSimplex2(slice<pt*, pt*> P) {

  // Maximize triangle area based on fixed xMin and xMax
  auto xx = minmax_element(P, [&](pt i, pt j) {return i[0]<j[0];});
  size_t xMin = xx.first - &P[0]; size_t xMax = xx.second - &P[0];
  pt x1 = P[xMin];
  pt x2 = P[xMax];
  auto y = max_element(P, [&](pt i, pt j) {
			    return crossProduct3d(x1-i, x2-i).length() <
			      crossProduct3d(x1-j, x2-j).length();
			  });
  size_t yApex = y - &P[0];
  pt y1 = P[yApex];

  // Maximize simplex volume
  pt area = crossProduct3d(x1-y1, x2-y1);
  auto z = max_element(P, [&](pt i, pt j) {
			    return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
			  });
  size_t zApex = z - &P[0];

  size_t c1 = xMin;
  size_t c2 = xMax;
  size_t c3 = yApex;
  size_t c4 = zApex;

  auto context = new _hull<linkedFacet3d, vertex3d>(P, c1, c2, c3, c4);

  return new conflictList(context);
}

template<class pt>
conflictList<linkedFacet3d, vertex3d> *makeInitialSimplex3(slice<pt*, pt*> P) {

  // Maximize triangle area based on fixed xMin and xMax
  size_t X[6]; // extrema
  auto xx = minmax_element(P, [&](pt i, pt j) {return i[0]<j[0];});
  X[0] = xx.first - &P[0]; X[1] = xx.second - &P[0];
  auto yy = minmax_element(P, [&](pt i, pt j) {return i[1]<j[1];});
  X[2] = yy.first - &P[0]; X[3] = yy.second - &P[0];
  auto zz = minmax_element(P, [&](pt i, pt j) {return i[2]<j[2];});
  X[4] = zz.first - &P[0]; X[5] = zz.second - &P[0];

  size_t xMin, xMax;
  if (P[X[1]][0]-P[X[0]][0] > P[X[3]][1]-P[X[2]][1] && P[X[1]][0]-P[X[0]][0] > P[X[5]][2]-P[X[4]][2]) {
    xMin = X[0]; xMax = X[1];
  } else if (P[X[3]][1]-P[X[2]][1] > P[X[1]][0]-P[X[0]][0] && P[X[3]][1]-P[X[2]][1] > P[X[5]][2]-P[X[4]][2]) {
    xMin = X[2]; xMax = X[3];
  } else {
    xMin = X[4]; xMax = X[5];
  }

  pt x1 = P[xMin];
  pt x2 = P[xMax];

  auto y = max_element(P, [&](pt i, pt j) {
			    return crossProduct3d(x1-i, x2-i).length() <
			      crossProduct3d(x1-j, x2-j).length();
			  });
  size_t yApex = y - &P[0];
  pt y1 = P[yApex];

  // Maximize simplex volume
  pt area = crossProduct3d(x1-y1, x2-y1);
  auto z = max_element(P, [&](pt i, pt j) {
			    return abs((y1-i).dot(area)) < abs((y1-j).dot(area));
			  });
  size_t zApex = z - &P[0];

  size_t c1 = xMin;
  size_t c2 = xMax;
  size_t c3 = yApex;
  size_t c4 = zApex;

  auto context = new _hull<linkedFacet3d, vertex3d>(P, c1, c2, c3, c4);

  return new conflictList(context);
}

////////////////////////////
// Main
////////////////////////////

template<class pt>
sequence<facet3d<pt>> incrementHull3dSerial(slice<pt*, pt*> P) {
  using facet3d = facet3d<pt>;

#ifdef WRITE
  {
    ofstream myfile;
    myfile.open("hull.txt", std::ofstream::trunc); myfile.close();
    myfile.open("point.txt", std::ofstream::trunc);
    for(size_t p=0; p<P.size(); ++p)
      myfile << P[p] << p << endl;
    myfile.close();
  }
#endif

  timer t; t.start();

  auto cg = makeInitialSimplex3(P);
  _hull<linkedFacet3d, vertex3d> *context = cg->context;

  double initTime = t.stop();

  size_t errors = 0;
  size_t round = 0;
  double apexTime = 0;
  double frontierTime = 0;
  double createTime = 0;
  double splitTime = 0;

  while (true) {
    t.start();

  loopStart:
    // Serial

#ifdef WRITE
    context->writeHull();
#endif

    //vertex3d apex = cg->randomApex();
    vertex3d apex = cg->furthestApex();

#ifdef VERBOSE
    cout << ">>>>>>>>>" << endl;
    //context->printHull();
    cout << "apex chosen = " << apex.attribute.i << endl;
#endif
    apexTime += t.get_next();

    if (apex.isEmpty()) break;

    round ++;

    auto frontier = cg->computeFrontier(apex);
    auto frontierEdges = get<0>(frontier);
    auto facetsBeneath = get<1>(frontier);

    frontierTime += t.get_next();

    for(size_t i=0; i<frontierEdges->size(); ++i) {
      auto nv = frontierEdges->at((i+1)%frontierEdges->size());
      auto cv = frontierEdges->at(i);
      if (cv.b != nv.a) {
	apex.attribute.seeFacet->clear();
	errors ++;
	goto loopStart;
      }
    }

#ifdef VERBOSE
    cout << "visible = " << *apex.attribute.seeFacet << ": ";
    for (auto x: *apex.attribute.seeFacet->seeList) {
      cout << x.attribute.i << " ";
    }
    cout << endl;

    cout << "frontier = ";
    for(auto e: *frontierEdges) {
      cout << e.a.attribute.i << "," << e.b.attribute.i << " ";
    }
    cout << endl;

    cout << "to delete = ";
    for(auto f: *facetsBeneath)
      cout << *f << " ";
    cout << endl;
#endif

    // Create new facets
    auto newFacets = sequence<linkedFacet3d*>(frontierEdges->size());

    for (size_t i=0; i<frontierEdges->size(); ++i) {
      _edge e = frontierEdges->at(i);
      newFacets[i] = makeFacet(e.a, e.b, apex);
    }

#ifdef VERBOSE
    cout << "to create = ";
    for (size_t i=0; i<frontierEdges->size(); ++i)
      cout << *(newFacets[i]) << " ";
    cout << endl;
#endif

    // Connect new facets
    for (size_t i=0; i<frontierEdges->size(); ++i) {
      linkFacet(newFacets[i],
		newFacets[(i+1)%frontierEdges->size()],
		frontierEdges->at(i).ff,
		newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		);
    }

    context->H = newFacets[0];

    createTime += t.get_next();

    cg->redistribute(facetsBeneath->cut(0, facetsBeneath->size()), make_slice(newFacets));

    splitTime += t.stop();

    // Delete existing facets
    for(int j=0; j<facetsBeneath->size(); ++j)
      delete facetsBeneath->at(j);

    delete frontierEdges;
    delete facetsBeneath;

  }

#ifndef SILENT
  cout << "init-time = " << initTime << endl;
  cout << "apex-time = " << apexTime << endl;
  cout << "frontier-time = " << frontierTime << endl;
  cout << "create-time = " << createTime << endl;
  cout << "split-time = " << splitTime << endl;
  cout << "#-rounds = " << round << endl;
  cout << "#-errors = " << errors << endl;
#endif

#ifdef VERBOSE
  cout << "hull-size = " << context->hullSizeDfs() << endl;
#endif

#ifdef WRITE
  context->writeHull();
#endif
  //free stuff

  auto result = sequence<facet3d>();
  context->getHull<pt>(result);
  return result;
}

template<class pt>
sequence<facet3d<pt>> incrementHull3d(slice<pt*, pt*> P, size_t numProc = 0) {
  using facet3d = facet3d<pt>;

  if (numProc == 0)
    numProc = parlay::num_workers();

#ifdef WRITE
  {
    ofstream myfile;
    myfile.open("hull.txt", std::ofstream::trunc); myfile.close();
    myfile.open("point.txt", std::ofstream::trunc);
    for(size_t p=0; p<P.size(); ++p)
      myfile << P[p] << p << endl;
    myfile.close();
  }
#endif

  timer t; t.start();

  auto cg = makeInitialSimplex3(P);

  _hull<linkedFacet3d, vertex3d> *context = cg->context;

  double initTime = t.stop();

  size_t round = 0;
  size_t serRound = 0;
  timer tr;
  double serTime = 0;
  double parTime = 0;

  double apexTime = 0;
  double frontierTime = 0;
  double splitTime = 0;

  while (true) {
    round ++;

    if (context->hullSize() < 64) {

      serRound ++;
      tr.start();
      t.start();

    loopStart:
      linkedFacet3d* f0 = context->hullSize() < 512 ? cg->facetWalk() : nullptr;
      vertex3d apex = cg->furthestApex(f0);

      apexTime += t.get_next();

      if (apex.isEmpty()) break;

      auto frontier = cg->computeFrontier(apex);
      auto frontierEdges = get<0>(frontier);
      auto facetsBeneath = get<1>(frontier);

      frontierTime += t.get_next();

      // Check for frontier error, usually caused by numerical errors in rare cases
      for(size_t i=0; i<frontierEdges->size(); ++i) {
	auto nv = frontierEdges->at((i+1)%frontierEdges->size());
	auto cv = frontierEdges->at(i);
	if (cv.b != nv.a) {
	  apex.attribute.seeFacet->clear();
	  goto loopStart;
	}
      }

      // Create new facets
      auto newFacets = sequence<linkedFacet3d*>(frontierEdges->size());

      for (size_t i=0; i<frontierEdges->size(); ++i) {
	_edge e = frontierEdges->at(i);
	newFacets[i] = makeFacet(e.a, e.b, apex);
      }

      // Connect new facets
      for (size_t i=0; i<frontierEdges->size(); ++i) {
	linkFacet(newFacets[i],
		  newFacets[(i+1)%frontierEdges->size()],
		  frontierEdges->at(i).ff,
		  newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		  );
      }

      context->H = newFacets[0];

      cg->redistribute(facetsBeneath->cut(0, facetsBeneath->size()), make_slice(newFacets));

      splitTime += t.stop();

      // Delete existing facets
      for(int j=0; j<facetsBeneath->size(); ++j)
	delete facetsBeneath->at(j);

      delete frontierEdges;
      delete facetsBeneath;

      serTime += tr.stop();
    } else {
      tr.start();
      t.start();

      // todo points chosen are close, which is bad
      sequence<vertex3d> apexes0 = cg->furthestApexes(numProc); // todo tune
      //sequence<vertex3d> apexes0 = cg->randomApexes(context->hullSize()); // todo tune
#ifdef VERBOSE
      size_t np = apexes0.size();
#endif

      apexTime += t.get_next();

      if (apexes0.size() <= 0) break;

      size_t numApex0 = apexes0.size();
      sequence<sequence<_edge<linkedFacet3d, vertex3d>>*> FE0(numApex0);
      sequence<sequence<linkedFacet3d*>*> FB0(numApex0);

      // Compute frontier and reserve facets
      parallel_for(0, numApex0, [&](size_t a) {
				 auto frontier = cg->computeFrontierAndReserve(apexes0[a]);

				 sequence<_edge<linkedFacet3d, vertex3d>>* frontierEdges = get<0>(frontier);
				 sequence<linkedFacet3d*>* facetsBeneath = get<1>(frontier);
				 FE0[a] = frontierEdges;
				 FB0[a] = facetsBeneath;
			       });

      frontierTime += t.get_next();

      // Check reservation for success, then reset reservation
      // Also check for numerical errors, set to unsuccessful if detected
      sequence<bool> success(numApex0);
      parallel_for(0, numApex0, [&](size_t a) {

				 // Check reservation
				 if (!cg->confirmReservation( apexes0[a], FB0[a]->cut(0, FB0[a]->size()) )) {
				   success[a] = false;
				 } else {
				   success[a] = true;
				 }

				 // Check for numerical errors in the frontier
				 for(size_t i=0; i<FE0[a]->size(); ++i) {
				   auto nv = FE0[a]->at((i+1)%FE0[a]->size());
				   auto cv = FE0[a]->at(i);
				   if (cv.b != nv.a) {
				     apexes0[a].attribute.seeFacet->clear();
				     success[a] = false;
				   }
				 }

				});
      parallel_for(0, numApex0, [&](size_t a) {
				 cg->resetReservation(apexes0[a], FB0[a]->cut(0, FB0[a]->size()));
			       });

      auto apexes = parlay::pack(make_slice(apexes0), success);
      auto FE = parlay::pack(make_slice(FE0), success);
      auto FB = parlay::pack(make_slice(FB0), success);
      size_t numApex = apexes.size();
#ifdef VERBOSE
      cout << numApex << "/" << np << endl;
#endif
      sequence<int> increase(numApex, 0);

      // Process the successful points
      parallel_for(0, numApex, [&](size_t a) {
				   // Create new facets
				   auto newFacets = sequence<linkedFacet3d*>(FE[a]->size());

				   for (size_t i=0; i<FE[a]->size(); ++i) {
				     _edge e = FE[a]->at(i);
				     newFacets[i] = makeFacet(e.a, e.b, apexes[a]);
				   }

				   // Connect new facets
				   for (size_t i=0; i<FE[a]->size(); ++i) {
				     linkFacet(newFacets[i],
					       newFacets[(i+1)%FE[a]->size()],
					       FE[a]->at(i).ff,
					       newFacets[(i-1+FE[a]->size())%FE[a]->size()]
					       );
				   }

				   context->H = newFacets[0]; // todo data race

				   cg->redistribute(FB[a]->cut(0, FB[a]->size()), make_slice(newFacets));

				   increase[a] = FE[a]->size() - FB[a]->size();
			       });

      //todo remove
      // if (!cg->checkReset()) {// checking this is expensive O(nh)
      //   cout << "not reset " << endl; abort();
      // }

      splitTime += t.stop();

      // Clean up
      parallel_for(0, numApex0, [&](size_t a) {
				 if (success[a]) {
				   // Delete existing facets
				   for(int j=0; j<FB0[a]->size(); ++j)
				     delete FB0[a]->at(j);
				 }

				 delete FE0[a];
				 delete FB0[a];
			       });

      parTime += tr.stop();
    } // End serial/parallel if

  }
#ifndef SILENT
  cout << "init-time = " << initTime << endl;
  cout << "apex-time = " << apexTime << endl;
  cout << "frontier-time = " << frontierTime << endl;
  cout << "create-split-time = " << splitTime << endl;
  cout << "#-rounds = " << round << endl;
  cout << "#-ser-rounds = " << serRound << endl;
  cout << "ser-time = " << serTime << endl;
  cout << "#-par-rounds = " << round - serRound << endl;
  cout << "par-time = " << parTime << endl;
  cout << "hull-size = " << context->hullSize() << endl;
#endif

#ifdef WRITE
  context->writeHull();
#endif

  auto result = sequence<facet3d>();
  context->getHull<pt>(result);
  return result;
}
