#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"
#include "euclideanMst/kruskal.h"
#include "euclideanMst/euclideanMst.h"

using namespace parlay;
using namespace pargeo;

template<int dim>
parlay::sequence<pargeo::wghEdge> pargeo::euclideanMst(parlay::sequence<pargeo::point<dim>> &S) {
  using nodeT = kdNode<dim, point<dim>>;
  using floatT = double;

  timer t;

  //nodeT* tree = buildKdt<dim, point<dim>>(S, true, true);
  nodeT* tree = buildKdTree<dim, point<dim>>(S, true, 1);

  cout << "build-time = " << t.get_next() << endl;

  auto pairs = wspdParallel(tree, 2);

  cout << "decomposition-time = " << t.get_next() << endl;

  cout << "pairs = " << pairs.size() << endl;

  struct wEdge {
    size_t u,v;
    floatT weight;
  };

  auto base = S.data();
  sequence<wEdge> edges = tabulate(pairs.size(), [&](size_t i) {
				  auto bcp = bccp<nodeT>(pairs[i].u, pairs[i].v);
				  wEdge e;
				  e.u = get<0>(bcp) - base;
				  e.v = get<1>(bcp) - base;
				  e.weight = get<2>(bcp);
				  return e;
				});

  cout << "bccp-time = " << t.get_next() << endl;

  auto edgeIds = kruskal(edges, S.size());

  cout << "kruskal-time = " << t.get_next() << endl;

  auto mstEdges = tabulate(edgeIds.size(),
			   [&](size_t i){
			     auto e = edges[edgeIds[i]];
			     return pargeo::wghEdge(e.u, e.v, e.weight);
			   });

  cout << "edges = " << mstEdges.size() << endl;
  return mstEdges;
}

template sequence<wghEdge> pargeo::euclideanMst<2>(sequence<point<2>> &);
template sequence<wghEdge> pargeo::euclideanMst<3>(sequence<point<3>> &);
template sequence<wghEdge> pargeo::euclideanMst<4>(sequence<point<4>> &);
template sequence<wghEdge> pargeo::euclideanMst<5>(sequence<point<5>> &);
template sequence<wghEdge> pargeo::euclideanMst<6>(sequence<point<6>> &);
template sequence<wghEdge> pargeo::euclideanMst<7>(sequence<point<7>> &);
