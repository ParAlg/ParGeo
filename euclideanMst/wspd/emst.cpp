#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/wspd.h"
#include "pargeo/kdTree.h"
#include "bccp.h"
#include "euclideanMst/euclideanMst.h"

using namespace parlay;
using namespace pargeo;

template<int dim>
parlay::sequence<pargeo::edge> pargeo::euclideanMst(parlay::sequence<pargeo::point<dim>> &S) {
  using nodeT = kdNode<dim, point<dim>>;
  using floatT = double;

  timer t;

  nodeT* tree = buildKdt<dim, point<dim>>(S, true, true);

  cout << "build-time = " << t.get_next() << endl;

  auto pairs = wspdParallel(tree, 2);

  cout << "decomposition-time = " << t.get_next() << endl;

  struct wEdge {
    size_t u,v;
    floatT w;
  };

  auto base = S.data();
  sequence<wEdge> edges = tabulate(pairs.size(), [&](size_t i) {
				  auto bcp = bccp<nodeT>(pairs[i].u, pairs[i].v);
				  wEdge e;
				  e.u = get<0>(bcp) - base;
				  e.v = get<1>(bcp) - base;
				  e.w = get<2>(bcp);
				  return e;
				});

  cout << "bccp-time = " << t.get_next() << endl;

  return sequence<pargeo::edge>();
}

template sequence<edge> pargeo::euclideanMst<2>(sequence<point<2>> &);
template sequence<edge> pargeo::euclideanMst<3>(sequence<point<3>> &);
template sequence<edge> pargeo::euclideanMst<4>(sequence<point<4>> &);
template sequence<edge> pargeo::euclideanMst<5>(sequence<point<5>> &);
template sequence<edge> pargeo::euclideanMst<6>(sequence<point<6>> &);
template sequence<edge> pargeo::euclideanMst<7>(sequence<point<7>> &);
