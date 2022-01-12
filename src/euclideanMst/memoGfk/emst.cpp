#include <tuple>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"
#include "euclideanMst/kruskal.h"
#include "euclideanMst/euclideanMst.h"
#include "wspdFilter.h"
#include "mark.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::emstInternal;

template<int dim>
parlay::sequence<pargeo::wghEdge> pargeo::euclideanMst(parlay::sequence<pargeo::point<dim>> &S) {
  using pointT = point<dim>;
  using nodeT = kdNode<dim, point<dim>>;
  using floatT = typename pointT::floatT;
  using pairT = wsp<nodeT>;
  using bcpT = tuple<pointT*, pointT*, floatT>;

  if (S.size() < 2) {
    throw std::runtime_error("need more than 2 points");
  }

  // timer t0;
  // t0.start();
  bool paraTree = true;

  //nodeT* tree = buildKdt<dim, point<dim>>(S, true, true);
  nodeT* tree = buildKdTree<dim, point<dim>>(S, true, 1);

  // cout << "build-tree-time = " << t0.get_next() << endl;

  floatT rhoLo = -0.1;
  floatT beta = 2;
  size_t numEdges = 0;

  floatT wspdTime = 0;
  floatT kruskalTime = 0;
  floatT markTime = 0;
  edgeUnionFind<long> UF(S.size());

  // t0.stop();

  while (UF.numEdge() < S.size() - 1) {

    // t0.start();

    floatT rhoHi;
    auto bccps = filterWspdParallel<nodeT>(beta, rhoLo, rhoHi, tree, &UF);

    // wspdTime += t0.get_next();

    // cout << "---" << endl;
    // cout << " beta = " << beta << endl;
    // cout << " rho = " << rhoLo << " -- " << rhoHi << endl;

    numEdges += bccps.size();

    if (bccps.size() <= 0) {
      beta *= 2;
      rhoLo = rhoHi;
      continue;}

    // cout << " edges = " << bccps.size() << endl;

    struct wEdge {
      size_t u,v;
      floatT weight;
    };

    auto base = S.data();
    sequence<wEdge> edges = tabulate(bccps.size(), [&](size_t i) {
	auto bcp = bccps[i];
	wEdge e;
	e.u = get<0>(bcp) - base;
	e.v = get<1>(bcp) - base;
	e.weight = get<2>(bcp);
	return e;
      });

    batchKruskal(edges, S.size(), UF);
    // cout << " mst-edges = " << UF.numEdge() << endl;
    // kruskalTime += t0.get_next();

    mark<nodeT, pointT, edgeUnionFind<long>>(tree, &UF, S.data());
    // markTime += t0.stop();

    beta *= 2;
    rhoLo = rhoHi;
  }

  // floatT sum = 0;
  // auto E = UF.getEdge();
  // for (auto e: E)
  //   sum += S[e.u].dist(S[e.v]);
  // cout << "edge-sum = " << sum << endl;

  // cout << "wspd-time = " << wspdTime << endl;
  // cout << "kruskal-time = " << kruskalTime << endl;
  // cout << "mark-time = " << markTime << endl;
  return UF.getEdge();
}

template sequence<wghEdge> pargeo::euclideanMst<2>(sequence<point<2>> &);
template sequence<wghEdge> pargeo::euclideanMst<3>(sequence<point<3>> &);
template sequence<wghEdge> pargeo::euclideanMst<4>(sequence<point<4>> &);
template sequence<wghEdge> pargeo::euclideanMst<5>(sequence<point<5>> &);
template sequence<wghEdge> pargeo::euclideanMst<6>(sequence<point<6>> &);
template sequence<wghEdge> pargeo::euclideanMst<7>(sequence<point<7>> &);
