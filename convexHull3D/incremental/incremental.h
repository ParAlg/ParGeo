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

#include <assert.h>
#include <stack>
#include <tuple>
#ifdef WRITE
#include <iostream>
#include <fstream>
#endif
#include "parlay/sequence.h"
#include "geometry/point.h"
#include "common/algebra.h"
#include "common/get_time.h"
#include "hull.h"
#include "split.h"

using namespace parlay;

////////////////////////////
// Primitives
////////////////////////////

/* Signed volume of an oriented tetrahedron (example below is positive).
     d
     o
    /|\
   / | o b
  o--o/
   a   c
*/
template <class pt>
inline typename pt::floatT signedVolume(pt a, pt b, pt c, pt d) {
  return (a-d).dot(crossProduct3d(b-a, c-a))/6;
}

////////////////////////////
// Facet
////////////////////////////

template <class pt>
struct facetAtt {
  _facet3d<pt, facetAtt<pt>> *abFacet;
  _facet3d<pt, facetAtt<pt>> *bcFacet;
  _facet3d<pt, facetAtt<pt>> *caFacet;
  slice<size_t*, size_t*> seeList;
  facetAtt() {}
};

template<class pt>
using linkedFacet3d = _facet3d<pt, facetAtt<pt>>;

template<class pt>
static std::ostream& operator<<(std::ostream& os, const linkedFacet3d<pt> v) {
  os << "(" << v.a << "," << v.b << "," << v.c << ")";
  return os;
}

// a, b, c are indices of vertices, can be input in any order
template<class pt>
linkedFacet3d<pt>* makeFacet(size_t a, size_t b, size_t c,
			     parlay::slice<pt*, pt*> P) {
  auto f = new linkedFacet3d<pt>(a, b, c, P);
  return f;
}

// Link f with ab, bc, ca; the edge matching is automatic -- input in any order
template<class fc>
void linkFacet(fc* f, fc* ab, fc* bc, fc* ca) {
  fc* F[3]; F[0]=ab; F[1]=bc; F[2]=ca;
  auto findFacet = [&](size_t v1, size_t v2) {
		     for(int i=0; i<3; ++i) {
		       if ((F[i]->a==v1 && F[i]->b==v2) || (F[i]->b==v1 && F[i]->a==v2) ||
			   (F[i]->b==v1 && F[i]->c==v2) || (F[i]->c==v1 && F[i]->b==v2) ||
			   (F[i]->c==v1 && F[i]->a==v2) || (F[i]->a==v1 && F[i]->c==v2)) {
			 return F[i];
		       }
		     }
		     cout << "Error, facet linking failure, abort." << endl;
		     abort();
		   };

  auto linkEdge = [&](fc* f1, fc* f2, size_t v1, size_t v2) {
		    if ( (f2->a==v1 && f2->b==v2) || (f2->a==v2 && f2->b==v1) )
		      f2->attribute.abFacet = f1;
		    else if ( (f2->b==v1 && f2->c==v2) || (f2->b==v2 && f2->c==v1) )
		      f2->attribute.bcFacet = f1;
		    else if ( (f2->c==v1 && f2->a==v2) || (f2->c==v2 && f2->a==v1) )
		      f2->attribute.caFacet = f1;
		    else {
		      cout << "Error, edge linking failure, abort." << endl;
		      abort();
		    }
		  };

  f->attribute.abFacet = findFacet(f->a, f->b);
  linkEdge(f, f->attribute.abFacet, f->a, f->b);
  f->attribute.bcFacet = findFacet(f->b, f->c);
  linkEdge(f, f->attribute.bcFacet, f->b, f->c);
  f->attribute.caFacet = findFacet(f->c, f->a);
  linkEdge(f, f->attribute.caFacet, f->c, f->a);
}

template<class facetT, class vertexT>
bool visible(facetT* f, size_t p, sequence<vertexT>* P) {
  if (signedVolume(P->at(f->a), P->at(f->b), P->at(f->c), P->at(p)) > 1e-7) // numeric knob
    return true;
  else
    return false;
}

template<class facetT, class vertexT>
bool visible(facetT* f, size_t p, slice<vertexT*, vertexT* > Q) {
  if (signedVolume(Q[f->a], Q[f->b], Q[f->c], Q[p]) > 1e-7) // numeric knob
    return true;
  else
    return false;
}

template <class pt>
bool visible(pt a, pt b, pt c, pt d) {
  return signedVolume(a, b, c, d) > 1e-7; // numeric knob
}

////////////////////////////
// Vertex
////////////////////////////

template<class facetT>
struct vertexAtt;

template<class facetT>
using vertex3d = pargeo::_point<3, double, double, vertexAtt<facetT>>;

template<class facetT>
struct vertexAtt {
  facetT *seeFacet;
  vertexAtt() {}
};

template<class facetT>
static std::ostream& operator<<(std::ostream& os, const vertex3d<facetT> v) {
  for (int i=0; i<v.dim; ++i)
    os << v.x[i] << " ";
  return os;
}

////////////////////////////
// Context
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
  size_t a, b;
  facetT* ff;
  facetT* fb;

  _edge(size_t _a, size_t _b, facetT* _ff, facetT* _fb) {
    a = _a; b = _b; ff = _ff; fb = _fb;}

  friend bool operator==(_edge e1, _edge e2) {
    return e1.a==e2.a && e1.b==e2.b;}
};

/* A context that maintains key states of the incremental algorithm
 */
template <class facetT, class vertexT>
struct _context {
  facetT* H;
  sequence<vertexT>* Q;

  _context(facetT* _H, sequence<vertexT>* _Q) {
    H = _H;
    Q = _Q;
  }

  /* Depth-first hull traversal (no facet repeat)
  */
  template <class F, class G, class H>
  void dfsFacet(facetT* start, F& fVisit, G& fDo, H& fStop) {
    sequence<facetT*> V;
    auto mark = [&](facetT* f) {V.push_back(f);};
    auto visited = [&](facetT* f) {
		     for (auto g: V) {
		       if (f == g) return true;
		     }
		     return false;};

    stack<_edge<facetT, vertexT>> S;

    S.emplace(start->b, start->a, start, nullptr);

    while (S.size() > 0) {
      _edge e = S.top(); S.pop();
      if (!visited(e.ff) && fVisit(e)) {
	fDo(e);
	mark(e.ff);

	S.emplace(e.ff->a, e.ff->b, e.ff->attribute.abFacet, e.ff);
	S.emplace(e.ff->c, e.ff->a, e.ff->attribute.caFacet, e.ff);
	S.emplace(e.ff->b, e.ff->c, e.ff->attribute.bcFacet, e.ff);
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
	  S.emplace(e.ff->a, e.ff->b, e.ff->attribute.abFacet, e.ff);
	  S.emplace(e.ff->c, e.ff->a, e.ff->attribute.caFacet, e.ff);
	  S.emplace(e.ff->b, e.ff->c, e.ff->attribute.bcFacet, e.ff);
	} else if (e.ff->b == e.b) {
	  S.emplace(e.ff->b, e.ff->c, e.ff->attribute.bcFacet, e.ff);
	  S.emplace(e.ff->a, e.ff->b, e.ff->attribute.abFacet, e.ff);
	  S.emplace(e.ff->c, e.ff->a, e.ff->attribute.caFacet, e.ff);
	} else if (e.ff->c == e.b) {
	  S.emplace(e.ff->c, e.ff->a, e.ff->attribute.caFacet, e.ff);
	  S.emplace(e.ff->b, e.ff->c, e.ff->attribute.bcFacet, e.ff);
	  S.emplace(e.ff->a, e.ff->b, e.ff->attribute.abFacet, e.ff);
	}
      }
      if (fStop()) break;
    }
  }

  /* Choose a random outside vertex visible to some facet
   */
  size_t randomApex() {
    if (H->attribute.seeList.size() > 0) return H->attribute.seeList[0];
    size_t apex = -1; // not unsigned
    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->attribute.seeList.size() > 0)
		   apex = e.ff->attribute.seeList[0];};
    auto fStop = [&]() {
		   if (apex != -1) return true;
		   else return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    return apex;
  }

  /* Choose the furthest outside vertex visible to some facet
   */
  size_t furthestApex() {
    //todo parallelize
    auto find = [&](slice<size_t*, size_t*> list, facetT* f) {
		  auto m = signedVolume(Q->at(f->a), Q->at(f->b), Q->at(f->c), Q->at(0));
		  size_t i = 0;
		  for (auto x: list) {
		    auto m2 = signedVolume(Q->at(f->a), Q->at(f->b), Q->at(f->c), Q->at(x));
		    if (m2 > m) {
		      m = m2; i = x;
		    }
		  }
		  return i;
		};

    if (H->attribute.seeList.size() > 0) return find(H->attribute.seeList, H);
    size_t apex = -1; // note this is unsigned
    auto fVisit = [&](_edge<facetT, vertexT> e) {return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (e.ff->attribute.seeList.size() > 0)
		   apex = find(e.ff->attribute.seeList, e.ff);};
    auto fStop = [&]() {
		   if (apex != -1) return true;
		   else return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    return apex;
  }

  /* Compute a frontier of edges in the clockwise order
      and facets to delete
   */
  tuple<sequence<_edge<facetT, vertexT>>*, sequence<facetT*>*> computeFrontier(size_t apex) {
    facetT* fVisible = Q->at(apex).attribute.seeFacet;

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
		    if (e.fb == nullptr || visible(e.fb, apex, Q))
		      return true;
		    else
		      return false;
  		  };

    auto fDo = [&](_edge<facetT, vertexT> e) {
		 // Include the facet for deletion if visible
		 bool seeff = visible(e.ff, apex, Q);
		 if ((seeff || e.fb == nullptr) && !facetVisited(e.ff))
		   facets->push_back(e.ff);

		 if (e.fb == nullptr) return; // Stop for the starting facet

		 // Include an edge joining a visible and an invisible facet as frontier
		 bool seefb = visible(e.fb, apex, Q);
		 if (seefb && !seeff) {
  		   frontier->emplace_back(e.a, e.b, e.ff, e.fb);
		 }
  	       };
    auto fStop = [&](){ return false;};

    dfsEdge(Q->at(apex).attribute.seeFacet, fVisit, fDo, fStop);
    return make_tuple(frontier, facets);
  }

  size_t hullSize(facetT* start=nullptr) {
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
    size_t hs = hullSize();
    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 if (checker && hullSize(e.ff) != hs) {
		   cout << " ..." << endl;
		   cout << "Error, inconsistency detected, abort" << endl;
		   cout << "Erroneous hull size = " << hullSize(e.ff) << endl;
		   abort();
		 }
		 if (!checker) cout << *e.ff << ":" << e.ff->attribute.seeList.size() << " ";
		 // if (!checker) {
		 //   cout << *e.ff;
		 //   if (e.ff->attribute.seeList.size() > 0) {
		 //     cout << ":";
		 //     for (auto x: e.ff->attribute.seeList) cout << x << " ";
		 //   } else
		 //     cout << " ";
		 // }
	       };
    auto fStop = [&]() { return false;};

    if (!checker) cout << "Hull DFS (" << hs << ")  = ";
    if (start) dfsFacet(start, fVisit, fDo, fStop);
    else dfsFacet(H, fVisit, fDo, fStop);
    if (!checker) cout << endl;
  }

#ifdef WRITE
  void writeHull() {
    ofstream myfile;
    myfile.open("hull.txt", std::ofstream::trunc);
    auto fVisit = [&](_edge<facetT, vertexT> e) { return true;};
    auto fDo = [&](_edge<facetT, vertexT> e) {
		 myfile << Q->at(e.ff->a) << endl;
		 myfile << Q->at(e.ff->b) << endl;
		 myfile << Q->at(e.ff->c) << endl;
	       };
    auto fStop = [&]() { return false;};
    dfsFacet(H, fVisit, fDo, fStop);
    myfile.close();
  }
#endif
};

////////////////////////////
// Incremental algorithm subroutines
////////////////////////////

template<class pt>
_context<linkedFacet3d<pt>, vertex3d<linkedFacet3d<pt>>> makeInitialHull(slice<pt*, pt*> P) {
  using linkedFacet3d = linkedFacet3d<pt>;
  using vertex3d = vertex3d<linkedFacet3d>;

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

// #ifdef WRITE
//   ofstream myfile;
//   myfile.open("hull.txt", std::ofstream::trunc);
//   myfile << P[xMin] << xMin << endl;
//   myfile << P[xMax] << xMax << endl;
//   myfile << P[yMin] << yMin << endl;
//   myfile << P[yMax] << yMax << endl;
//   myfile << P[zMin] << zMin << endl;
//   myfile << P[zMax] << zMax << endl;
//   myfile.close();
// #endif

  // Make initial facets
  auto f0 = makeFacet(xMin, yMin, zMin, P);
  auto f1 = makeFacet(xMax, yMin, zMin, P);
  auto f2 = makeFacet(xMin, yMax, zMin, P);
  auto f3 = makeFacet(xMax, yMax, zMin, P);
  auto f4 = makeFacet(xMin, yMin, zMax, P);
  auto f5 = makeFacet(xMax, yMin, zMax, P);
  auto f6 = makeFacet(xMin, yMax, zMax, P);
  auto f7 = makeFacet(xMax, yMax, zMax, P);

  linkFacet(f0, f2, f1, f4);
  linkFacet(f1, f5, f0, f3);
  linkFacet(f2, f6, f3, f0);
  linkFacet(f3, f1, f2, f7);
  linkFacet(f4, f0, f5, f6);
  linkFacet(f5, f7, f4, f1);
  linkFacet(f6, f4, f7, f2);
  linkFacet(f7, f3, f6, f5);
  // Initialize points with visible facet link
  auto Q = new sequence<vertex3d>(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      for(int j=0; j<P[i].dim; ++j)
				Q->at(i)[j] = P[i][j];
			    });

  // Assign each point to one of the visible facets
  auto flag = sequence<int>(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      if (visible(P[f0->a], P[f0->b], P[f0->c], P[i])) {
				flag[i] = 0; Q->at(i).attribute.seeFacet = f0;
			      } else if (visible(P[f1->a], P[f1->b], P[f1->c], P[i])) {
				flag[i] = 1; Q->at(i).attribute.seeFacet = f1;
			      } else if (visible(P[f2->a], P[f2->b], P[f2->c], P[i])) {
				flag[i] = 2; Q->at(i).attribute.seeFacet = f2;
			      } else if (visible(P[f3->a], P[f3->b], P[f3->c], P[i])) {
				flag[i] = 3; Q->at(i).attribute.seeFacet = f3;
			      } else if (visible(P[f4->a], P[f4->b], P[f4->c], P[i])) {
				flag[i] = 4; Q->at(i).attribute.seeFacet = f4;
			      } else if (visible(P[f5->a], P[f5->b], P[f5->c], P[i])) {
				flag[i] = 5; Q->at(i).attribute.seeFacet = f5;
			      } else if (visible(P[f6->a], P[f6->b], P[f6->c], P[i])) {
				flag[i] = 6; Q->at(i).attribute.seeFacet = f6;
			      } else if (visible(P[f7->a], P[f7->b], P[f7->c], P[i])) {
				flag[i] = 7; Q->at(i).attribute.seeFacet = f7;
			      } else {
				flag[i] = 8; Q->at(i).attribute.seeFacet = nullptr;
			      }
			    });

  // Partition points based on the facets assigned to
  auto I = sequence<size_t>(P.size());
  auto I2 = new sequence<size_t>(P.size()); //todo free
  parallel_for(0, P.size(), [&](size_t i){I[i] = i;});
  auto splits = split_k(9, make_slice(I), make_slice(I2->begin(), I2->end()), flag);

  size_t n = scan_inplace(make_slice(splits), addm<size_t>());
  splits.push_back(n);

  // Assign points to the respecive facets
  f0->attribute.seeList = I2->cut(splits[0], splits[1]);
  f1->attribute.seeList = I2->cut(splits[1], splits[2]);
  f2->attribute.seeList = I2->cut(splits[2], splits[3]);
  f3->attribute.seeList = I2->cut(splits[3], splits[4]);
  f4->attribute.seeList = I2->cut(splits[4], splits[5]);
  f5->attribute.seeList = I2->cut(splits[5], splits[6]);
  f6->attribute.seeList = I2->cut(splits[6], splits[7]);
  f7->attribute.seeList = I2->cut(splits[7], splits[8]);

  return _context<linkedFacet3d, vertex3d>(f0, Q);
}

////////////////////////////
// Main
////////////////////////////

template<class pt>
sequence<facet3d<pt>> incrementHull3d(slice<pt*, pt*> P) {
  using facet3d = facet3d<pt>;
  using linkedFacet3d = linkedFacet3d<pt>;
  using vertex3d = vertex3d<linkedFacet3d>;

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

  _context<linkedFacet3d, vertex3d> context = makeInitialHull(P);

  cout << "init-time = " << t.get_next() << endl;

#ifdef VERBOSE
  size_t round = 0;
  size_t simpleSplit = 0;
  size_t allocSplit = 0;
#endif

  while (true) {

    // Serial

    //size_t apex = context.randomApex();
    size_t apex = context.furthestApex();
    if (apex == -1) break;

    round ++;

#ifdef VERBOSE
    // cout << ">>>>>>>>>> round" << endl;
    // context.printHull();
    // cout << "apex = " << apex << endl;
#endif
#ifdef WRITE
    context.writeHull();
#endif

    auto frontier = context.computeFrontier(apex);
    auto frontierEdges = get<0>(frontier);
    auto facetsBeneath = get<1>(frontier);

#ifdef VERBOSE
    // cout << "frontier = ";
    // for(auto e: *frontierEdges)
    //   cout << e.a << "," << e.b << " ";
    // cout << endl;

    // cout << "to delete = ";
    // for(auto f: *facetsBeneath)
    //   cout << *f << " ";
    // cout << endl;
#endif

    // Create new facets
    linkedFacet3d* newFacets[frontierEdges->size()];

    for (size_t i=0; i<frontierEdges->size(); ++i) {
      _edge e = frontierEdges->at(i);
      newFacets[i] = makeFacet(e.a, e.b, apex, P);
    }

#ifdef VERBOSE
    // cout << "to create = ";
    // for (size_t i=0; i<frontierEdges->size(); ++i)
    //   cout << *(newFacets[i]) << " ";
    // cout << endl;
#endif

    // Connect new facets
    for (size_t i=0; i<frontierEdges->size(); ++i) {
      linkFacet(newFacets[i],
		newFacets[(i+1)%frontierEdges->size()],
		frontierEdges->at(i).ff,
		newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		);
    }

    context.H = newFacets[0];

    // Redistribute the outside points
    if (facetsBeneath->size() == 1) {
      // Only one facet affected, simply distribute among three facets

#ifdef VERBOSE
      simpleSplit ++;
#endif

      linkedFacet3d* fdel = facetsBeneath->at(0);
      size_t fn = fdel->attribute.seeList.size();
      auto flag = sequence<int>(fn);
      parallel_for(0, fn, [&](size_t i) {
			    flag[i] = 3;
			    context.Q->at(fdel->attribute.seeList[i]).attribute.seeFacet = nullptr;
			    for (int j=0; j<3; ++j) {
			      if (visible(newFacets[j], fdel->attribute.seeList[i], P)) {
				flag[i] = j;
				context.Q->at(fdel->attribute.seeList[i]).attribute.seeFacet = newFacets[j];
				break;
			      }
			    }
			  });

      // Partition points based on the facets assigned to
      auto tmpBuffer = sequence<size_t>(fn);
      auto splits = split_k(4, fdel->attribute.seeList, make_slice(tmpBuffer), flag);

      size_t fn2 = scan_inplace(make_slice(splits), addm<size_t>());
      splits.push_back(fn2);

      parallel_for(0, fn, [&](size_t i) {fdel->attribute.seeList[i] = tmpBuffer[i];});

      for (int j=0; j<3; ++j)
	newFacets[j]->attribute.seeList =
	  fdel->attribute.seeList.cut(splits[j], splits[j+1]);

    } else {
      // Multiple facets affected, new memory needed

#ifdef VERBOSE
      allocSplit ++;
#endif

      int nf = facetsBeneath->size();
      int nnf = frontierEdges->size();

      size_t fn = 0;
      for(int j=0; j<nf; ++j) {
	fn += facetsBeneath->at(j)->attribute.seeList.size();
      }

      auto tmpBuffer = sequence<size_t>(fn);
      fn = 0;
      for(int j=0; j<nf; ++j) {
	parallel_for(0, facetsBeneath->at(j)->attribute.seeList.size(),
		     [&](size_t x){
		       tmpBuffer[fn+x] = facetsBeneath->at(j)->attribute.seeList[x];
		     });
	fn += facetsBeneath->at(j)->attribute.seeList.size();
      }

      auto flag = sequence<int>(fn);
      parallel_for(0, fn, [&](size_t i) {
			    flag[i] = nnf;
			    context.Q->at(tmpBuffer[i]).attribute.seeFacet = nullptr;
			    for (int j=0; j<nnf; ++j) {
			      if (visible(newFacets[j], tmpBuffer[i], P)) {
				flag[i] = j;
				context.Q->at(tmpBuffer[i]).attribute.seeFacet = newFacets[j];
				break;
			      }
			    }
			  });

      // Partition points based on the facets assigned to
      auto tmpBuffer2 = new sequence<size_t>(fn);//todo free
      auto splits = split_k(nnf+1, make_slice(tmpBuffer), tmpBuffer2->cut(0, tmpBuffer2->size()), flag);

      size_t fn2 = scan_inplace(make_slice(splits), addm<size_t>());
      splits.push_back(fn2);

      for (int j=0; j<nnf; ++j)
	newFacets[j]->attribute.seeList =
	  tmpBuffer2->cut(splits[j], splits[j+1]);
    }

    // Delete existing facets
    for(int j=0; j<facetsBeneath->size(); ++j)
      delete facetsBeneath->at(j);

  }
  cout << "increment-time = " << t.get_next() << endl;

#ifdef VERBOSE
  cout << "hull-size = " << context.hullSize() << endl;
  cout << "#rounds = " << round << endl;
  cout << "simple-splits = " << simpleSplit << endl;
  cout << "alloc-splits = " << allocSplit << endl;
#endif

#ifdef WRITE
  context.writeHull();
#endif

  // todo
  auto dummy = sequence<facet3d>(); dummy.reserve(1);
  return dummy;
}
