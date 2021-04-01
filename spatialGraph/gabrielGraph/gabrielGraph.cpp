#include <algorithm>
#include <tuple>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "geometry/point.h"
#include "common/geometry.h"
#include "common/get_time.h"
#include "spatialGraph.h"
#include "incremental/delaunay.h"


template<int dim>
parlay::sequence<edge> spatialGraph(parlay::sequence<pargeo::point<dim>> &P) {
  using namespace parlay;
  using namespace std;

  if (dim != 2)
    throw std::runtime_error("Error, gabriel graph only supports 2 dimensional inputs.");

  timer t; t.start();
  // A pbbsbench point data structure
  sequence<point2d<double>> P2d(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      P2d[i].x = P[i][0];
			      P2d[i].y = P[i][1];
			    });
  cout << "data-convert-time = " << t.get_next() << endl;

  triangles<point2d<double>> Tri = delaunay(P2d);
  cout << "triangulation-time = " << t.get_next() << endl;

  using pt = pargeo::point<dim>;
  using edgepair = tuple<edge,bool>;

  size_t nt = Tri.numTriangles();
  sequence<edgepair> edges(nt*3+1);
  edges[edges.size()-1] = make_tuple(edge(), false);

  // Process a triangle
  auto processTri = [&](tri &T, size_t idx) {
		      for(int i=0; i<3; ++i) {
			size_t u = T[(i-1)%3];
			size_t v = T[(i+1)%3];
			pt c = (P[u] + P[v])/2;
			typename pt::floatT rq = (P[u] - P[v]).length()/2;
			typename pt::floatT dis = (P[i] - c).length();
			if (dis <= rq) { // not ok
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

  // Keep valid edges
  sequence<size_t> flag(nt*3+1);
  parallel_for(0, edges.size(), [&](size_t i) {
				  if ( get<1>(edges[i]) ) {
				    if ( get<0>(edges[i])!=get<0>(edges[i+1]) ) {
				      flag[i] = true;
				    }
				    else flag[i] = false;
				  }
				  else flag[i] = false;
				});

  size_t ne = scan_inplace(flag.cut(0,nt*3));
  flag[flag.size()-1] = ne;

  sequence<edge> edges2(ne);
  parallel_for(0, nt*3, [&](size_t i) {
			if (flag[i] != flag[i+1]) {
			  edges2[flag[i]] = get<0>(edges[i]);
			}
  		      });

  // parallel_for(0, edges2.size(), [&](size_t i) {
  // 				  auto e = edges2[i];
  // 				  cout << "(" << e.u << "," << e.v << ")" << " ";
  // 				});
  // cout << endl;

  cout << "graph-gen-time = " << t.stop() << endl;

  cout << "#-triangles = " << nt << endl;
  cout << "#gabriel-edges = " << ne << endl;
  return edges2;
}

template parlay::sequence<edge> spatialGraph<2>(parlay::sequence<pargeo::point<2>> &);
template parlay::sequence<edge> spatialGraph<3>(parlay::sequence<pargeo::point<3>> &);
