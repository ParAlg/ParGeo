#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "pargeo/point.h"
#include "pargeo/kdTree.h"
#include "pargeo/wspd.h"
#include "pargeo/getTime.h"
#include "spatialGraph/spatialGraph.h"

template<int dim>
parlay::sequence<pargeo::edge> pargeo::spanner(parlay::sequence<pargeo::point<dim>> &S, double t) {
  using namespace std;
  using namespace pargeo;
  using namespace parlay;

  using pointT = point<dim>;
  using nodeT = kdNode<dim, pointT>;
  using pairT = wsp<nodeT>;

  // cout << t << "-spanner of " << S.size() << ", dim " << dim << " points" << endl;
  if (S.size() < 2) abort();

  double s = 4*((double)t+1)/((double)t-1);
  // cout << "separation-constant = " << s << endl;

  if (s > 10) std::cout << "WARNING: very large separation constant, try a smaller t if crashed\n";

  timer t0;
  t0.start();

  nodeT* tree = buildKdTree<dim, point<dim>>(S, true, 1);
  // cout << "build-tree-time = " << t0.get_next() << endl;

  auto wspd = wspdParallel<dim>(tree, s);
  // cout << "decomposition-time = " << t0.get_next() << endl;

  auto edges = sequence<edge>(wspd.size());
  pointT* base = S.data();
  parallel_for(0, wspd.size(), [&](size_t i) {
      size_t pt1 = wspd.at(i).u->getItem(0) - base;
      size_t pt2 = wspd.at(i).v->getItem(0) - base;
      edges[i] = edge(pt1, pt2);
    });

  // cout << "edge-count = " << edges.size() << endl;

  return edges;
}

template parlay::sequence<pargeo::edge> pargeo::spanner<2>(parlay::sequence<pargeo::point<2>> &, double);
template parlay::sequence<pargeo::edge> pargeo::spanner<3>(parlay::sequence<pargeo::point<3>> &, double);
template parlay::sequence<pargeo::edge> pargeo::spanner<4>(parlay::sequence<pargeo::point<4>> &, double);
template parlay::sequence<pargeo::edge> pargeo::spanner<5>(parlay::sequence<pargeo::point<5>> &, double);
template parlay::sequence<pargeo::edge> pargeo::spanner<6>(parlay::sequence<pargeo::point<6>> &, double);
template parlay::sequence<pargeo::edge> pargeo::spanner<7>(parlay::sequence<pargeo::point<7>> &, double);
