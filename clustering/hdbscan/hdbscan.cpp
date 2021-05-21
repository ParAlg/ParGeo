#include <tuple>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/wspd.h"
#include "pargeo/kdTree.h"
#include "pargeo/kdTreeKnn.h"
#include "pargeo/bccp.h"
#include "pargeo/kruskal.h"
#include "clustering/hdbscan.h"
#include "wspdFilter.h"
#include "mark.h"

using namespace std;
using namespace parlay;
using namespace pargeo;
using namespace pargeo::hdbscanInternal;

template<int dim>
parlay::sequence<pargeo::edge> pargeo::hdbscan(parlay::sequence<pargeo::point<dim>> &S, size_t minPts) {
  using pointT = point<dim>;
  using nodeT = kdNode<dim, point<dim>>;
  using floatT = typename pointT::floatT;
  using pairT = wsp<nodeT>;
  using bcpT = tuple<pointT*, pointT*, floatT>;

  if (S.size() < 2) {
    throw std::runtime_error("need more than 2 points");
  }

  timer t0;
  t0.start();

  nodeT* tree = buildKdt<dim, point<dim>>(S, true, true);

  cout << "build-tree-time = " << t0.get_next() << endl;

  // todo return distances
  sequence<size_t> nns = kdTreeKnn<dim, pointT>(S, minPts, tree, true);

  sequence<floatT> coreDist = sequence<floatT>(S.size());
  parallel_for (0, S.size(), [&](intT i) {
			       coreDist[i] = S[nns[i*minPts + minPts-1]].dist(S[i]);
			     });

  cout << "core-dist-time = " << t0.get_next() << endl;

  sequence<floatT> cdMin = sequence<floatT>(tree->size() * 2);
  sequence<floatT> cdMax = sequence<floatT>(tree->size() * 2);
  parallel_for(0, tree->size()*2, [&](intT i) {
      cdMin[i] = std::numeric_limits<floatT>::max();
      cdMax[i] = std::numeric_limits<floatT>::lowest();
    });
  hdbscanInternal::nodeCD(tree, coreDist, cdMin, cdMax, tree, S.data());

  floatT rhoLo = -0.1;
  floatT beta = 2;
  size_t numEdges = 0;

  floatT wspdTime = 0;
  floatT kruskalTime = 0;
  floatT markTime = 0;
  edgeUnionFind<long> UF(S.size());

  t0.stop();

  while (UF.numEdge() < S.size() - 1) {

    t0.start();

    floatT rhoHi;
    auto bccps = filterWspdParallel<nodeT>(beta, rhoLo, rhoHi, tree, &UF,
					   coreDist, cdMin, cdMax, S.data());

    wspdTime += t0.get_next();

    cout << "---" << endl;
    cout << " beta = " << beta << endl;
    cout << " rho = " << rhoLo << " -- " << rhoHi << endl;

    numEdges += bccps.size();

    if (bccps.size() <= 0) {
      beta *= 2;
      rhoLo = rhoHi;
      continue;}

    cout << " edges = " << bccps.size() << endl;

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
    cout << " mst-edges = " << UF.numEdge() << endl;
    kruskalTime += t0.get_next();

    mark<nodeT, pointT, edgeUnionFind<long>>(tree, &UF, S.data());
    markTime += t0.stop();

    beta *= 2;
    rhoLo = rhoHi;
  }

  // floatT sum = 0;
  // auto E = UF.getEdge();
  // for (auto e: E) {
  //   floatT w = S[e.u].dist(S[e.v]);
  //   w = max(w, coreDist[e.u]);
  //   w = max(w, coreDist[e.v]);
  //   sum += w;
  // }
  // cout << "edge-sum = " << sum << endl;

  cout << "wspd-time = " << wspdTime << endl;
  cout << "kruskal-time = " << kruskalTime << endl;
  cout << "mark-time = " << markTime << endl;
  return UF.getEdge();
}

template sequence<edge> pargeo::hdbscan<2>(sequence<point<2>> &, size_t);
template sequence<edge> pargeo::hdbscan<3>(sequence<point<3>> &, size_t);
template sequence<edge> pargeo::hdbscan<4>(sequence<point<4>> &, size_t);
template sequence<edge> pargeo::hdbscan<5>(sequence<point<5>> &, size_t);
template sequence<edge> pargeo::hdbscan<6>(sequence<point<6>> &, size_t);
template sequence<edge> pargeo::hdbscan<7>(sequence<point<7>> &, size_t);
