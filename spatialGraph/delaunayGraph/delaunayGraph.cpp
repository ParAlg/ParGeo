#include <algorithm>
#include <tuple>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "spatialGraph/spatialGraph.h"
#include "delaunayTriangulation/delaunay.h"
#include "delaunayTriangulation/geometry.h"

template<int dim>
parlay::sequence<pargeo::edge> pargeo::delaunayGraph(parlay::sequence<pargeo::point<dim>> &P) {
  using namespace pbbsbench;
  using namespace parlay;
  using namespace std;

  if (dim != 2)
    throw std::runtime_error("Error, delaunay graph only supports 2 dimensional inputs.");
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

  size_t n = P.size();
  size_t nt = Tri.numTriangles();
  sequence<edge> edges(nt*3+1);
  edges[nt*3] = edge();

  // Process a triangle
  auto processTri = [&](tri &T, size_t idx) {
		      for(int i=0; i<3; ++i) {
			size_t u = T[(i+1)%3];
			size_t v = T[(i+2)%3];

			// Keep non-boundary edges
			if (u < n && v < n)
			  edges[idx*3 + i] = edge(u,v);
			else
			  edges[idx*3 + i] = edge();

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
				  if ( edges[i] == edges[i+1] && !edges[i].isEmpty() )
				    flag[i] = 1;
				  else
				    flag[i] = 0;
				});

  sequence<edge> edges2 = pack(edges, flag);

#ifndef SILENT
  cout << "graph-gen-time = " << t.stop() << endl;
  cout << "#-triangles = " << nt << endl;
  cout << "#delaunay-edges = " << edges2.size() << endl;
#endif
  return edges2;
}

template parlay::sequence<pargeo::edge> pargeo::delaunayGraph<2>(parlay::sequence<pargeo::point<2>> &);
template parlay::sequence<pargeo::edge> pargeo::delaunayGraph<3>(parlay::sequence<pargeo::point<3>> &);
