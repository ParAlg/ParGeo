#include <algorithm>
#include <tuple>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "geometry/point.h"
#include "common/geometry.h"
#include "common/get_time.h"
#include "spatialGraph.h"
#include "delaunay.h"

template<int dim>
parlay::sequence<pargeo::edge> pargeo::gabrielGraph(parlay::sequence<pargeo::point<dim>> &P) {
  using namespace parlay;
  using namespace std;
  using namespace pbbsbench;

  if (dim != 2)
    throw std::runtime_error("Error, gabriel graph only supports 2 dimensional inputs.");
#ifndef SILENT
  timer t; t.start();
#endif
  // A pbbsbench point data structure
  sequence<point2d<double>> P2d(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      P2d[i].x = P[i][0];
			      P2d[i].y = P[i][1];
			    });
#ifndef SILENT
  cout << "data-convert-time = " << t.get_next() << endl;
#endif
  triangles<point2d<double>> Tri = delaunay(P2d);
#ifndef SILENT
  cout << "triangulation-time = " << t.get_next() << endl;
#endif
  using pt = pargeo::point<dim>;
  using edgepair = tuple<edge,bool>;

  size_t n = P.size();
  size_t nt = Tri.numTriangles();
  sequence<edgepair> edges(nt*3+1);
  edges[edges.size()-1] = make_tuple(edge(), false);

  // Process a triangle
  auto processTri = [&](tri &T, size_t idx) {
		      for(int i=0; i<3; ++i) {
			// Testing if x violates u,v
			size_t x = T[i];
			size_t u = T[(i+1)%3];
			size_t v = T[(i+2)%3];

			// Discard irrelevant edges of boundary triangle edges
			if (u >= n || v >= n) {
			  edges[idx*3 + i] = make_tuple(edge(u,v), false);
			  continue;
			}

			pt c = (P[u] + P[v])/2;
			typename pt::floatT rq = (P[u] - P[v]).length()/2;
			typename pt::floatT dis = (P[x] - c).length();
			if (x < n && dis <= rq) { // not gabriel
			  edges[idx*3 + i] = make_tuple(edge(u,v), false);
			} else { // ok
			  edges[idx*3 + i] = make_tuple(edge(u,v), true);
			}
		      }
		    };

  // Process triangles
  parallel_for(0, nt, [&](size_t i) {
			processTri(Tri.T[i], i);
		      });

  // Group edges
  sort_inplace(edges, [&](edgepair e1, edgepair e2) {
		return get<0>(e1).u == get<0>(e2).u ?
		  (get<0>(e1).v < get<0>(e2).v):
		  (get<0>(e1).u < get<0>(e2).u);
	      });

  // parallel_for(0, edges.size(), [&](size_t i) {
  // 				  auto e = get<0>(edges[i]);
  // 				  cout << "(" << e.u << "," << e.v << ")" << get<1>(edges[i]) << " ";
  // 				});
  // cout << endl;

  // Keep valid edges
  //  - appears twice && both are valid in their own triangles
  sequence<size_t> flag(nt*3+1);
  parallel_for(0, edges.size(), [&](size_t i) {
				  if ( get<0>(edges[i]) == get<0>(edges[i+1]) &&
				       get<1>(edges[i]) &&
				       get<1>(edges[i+1]) ) {
				    flag[i] = 1;
				  } else {
				    flag[i] = 0;
				  }
				});

  size_t ne = scan_inplace(flag.cut(0,nt*3));
  flag[nt*3] = ne;

  sequence<pargeo::edge> edges2(ne);
  parallel_for(0, nt*3, [&](size_t i) {
			if (flag[i] != flag[i+1]) {
			  edges2[flag[i]] = get<0>(edges[i]);
			}
  		      });

#ifndef SILENT
  cout << "graph-gen-time = " << t.stop() << endl;
  cout << "#-triangles = " << nt << endl;
  cout << "#gabriel-edges = " << ne << endl;
#endif
  return edges2;
}

template parlay::sequence<pargeo::edge> pargeo::gabrielGraph<2>(parlay::sequence<pargeo::point<2>> &);
