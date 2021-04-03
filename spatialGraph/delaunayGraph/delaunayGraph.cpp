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
    throw std::runtime_error("Error, delaunay graph only supports 2 dimensional inputs.");

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

  size_t n = P.size();
  size_t nt = Tri.numTriangles();
  sequence<edge> edges(nt*3+1);
  edges[nt*3] = edge();

  // Process a triangle
  auto processTri = [&](tri &T, size_t idx) {
		      for(int i=0; i<3; ++i) {
			size_t x = T[i];
			size_t u = T[(i+1)%3];
			size_t v = T[(i+2)%3];

			// Keep non-boundary edges
			if (u < n && v < n)
			  edges[idx*3 + i] = edge(u,v);

		      }
		    };

  // Process triangles
  parallel_for(0, nt, [&](size_t i) {
			processTri(Tri.T[i], i);
		      });

  // Group edges
  sort_inplace(edges, [&](edge e1, edge e2) {
			return e1.u == e2.u ?
			  (e1.v < e2.v) : (e1.u < e2.u);
		      });

  // Only keep each edge once
  sequence<size_t> flag(nt*3);
  parallel_for(0, edges.size(), [&](size_t i) {
				  if ( edges[i] == edges[i+1] )
				    flag[i] = 1;
				  else
				    flag[i] = 0;
				});

  sequence<edge> edges2 = pack(edges, flag);

  cout << "graph-gen-time = " << t.stop() << endl;

  cout << "#-triangles = " << nt << endl;
  cout << "#delaunay-edges = " << edges2.size() << endl;
  return edges2;
}

template parlay::sequence<edge> spatialGraph<2>(parlay::sequence<pargeo::point<2>> &);
template parlay::sequence<edge> spatialGraph<3>(parlay::sequence<pargeo::point<3>> &);
