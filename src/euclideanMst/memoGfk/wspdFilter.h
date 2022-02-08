// This code is part of the project "Fast Parallel Algorithms for Euclidean
// Minimum Spanning Tree and Hierarchical Spatial Clustering"
// Copyright (c) 2021 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
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

#include <atomic>
#include <tuple>
#include <limits>
#include "kdTree/kdTree.h"
#include "pargeo/parBuffer.h"
#include "pargeo/getTime.h"
#include "parlay/utilities.h"

namespace pargeo {
namespace emstInternal {

template<class nodeT, class UF>
struct rhoUpdateSerial {
  using floatT = double;

  floatT rho;
  floatT beta;
  UF *uf;

  rhoUpdateSerial(floatT betaa, UF *uff) :
    beta(betaa), uf(uff) {
    rho = std::numeric_limits<floatT>::max();}

  void run(nodeT *u, nodeT *v) {
    floatT myDist = nodeDistance(u, v);
    if (myDist < rho) rho = myDist;
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) return false; // filtering todo
    if (u->size() + v->size() <= beta) return false; // not in E_u, not considered for rho
    if (nodeDistance(u, v) >= rho) return false; // no subsequent finds can update rho
    return true;
  }

  bool start(nodeT *u) {
    if (u->size() > beta) {
      return true;
    } else {
      return false;// if node size < beta, so would children
    }
  }

  floatT getRho() { return rho;}

  bool wellSeparated(nodeT* u, nodeT* v, floatT s) {return geomWellSeparated(u, v, s);}
};

template<class nodeT, class UF>
struct wspGetSerial {
  using floatT = double;
  using bcpT = std::tuple<typename nodeT::objT*,
			  typename nodeT::objT*,
			  typename nodeT::objT::floatT>;

  floatT rhoLo;
  floatT rhoHi;
  floatT beta;
  parlay::sequence<bcpT> out;
  nodeT *tree;
  UF *uf;

  wspGetSerial(floatT betaa, floatT rhoLoo, floatT rhoHii, nodeT *treee, UF *uff) :
    beta(betaa), rhoLo(rhoLoo), rhoHi(rhoHii), tree(treee), uf(uff) {
    out = parlay::sequence<bcpT>();
  }

  void run(nodeT *u, nodeT *v) {
    auto bcp = bichromaticClosestPair(u, v);
    if (u->size() + v->size() <= beta &&
        std::get<2>(bcp) >= rhoLo &&
        std::get<2>(bcp) < rhoHi) {
      out.push_back(bcp);
    }
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) {return false;} // todo need to add id to tree nodes
    if (nodeDistance(u, v) >= rhoHi) return false; // too separated to be considered
    if (nodeFarDistance(u, v) < rhoLo) return false; // too close to be considered
    return true;
  }

  bool start(nodeT *u) {
    if (u->diag() >= rhoLo) {
      return true;
    } else {
      return false;
    }
  }

  bool wellSeparated(nodeT* u, nodeT* v, floatT s) {return geomWellSeparated(u, v, s);}
  parlay::sequence<bcpT> collect() { return out; };
};

template <class nodeT, class UF>
parlay::sequence<std::tuple<typename nodeT::objT*,
			    typename nodeT::objT*,
			    typename nodeT::objT::floatT>>
filterWspdSerial(double t_beta,
		 double t_rho_lo,
		 double& t_rho_hi,
		 nodeT *t_kdTree,
		 UF *t_mst) {
  using floatT = double;
  using objT = typename nodeT::objT;
  using bcpT = std::tuple<objT*, objT*, floatT>;

  auto myRho = rhoUpdateSerial<nodeT, UF>(t_beta, t_mst);

  pargeo::kdTree::computeWspdSerial<nodeT, rhoUpdateSerial<nodeT, UF>>(t_kdTree, &myRho);

  auto mySplitter = wspGetSerial<nodeT, UF>(t_beta, t_rho_lo, myRho.getRho(), t_kdTree, t_mst);

  pargeo::kdTree::computeWspdSerial<nodeT, wspGetSerial<nodeT, UF>>(t_kdTree, &mySplitter);

  t_rho_hi = myRho.getRho();
  return mySplitter.collect();
}

template<class nodeT, class UF>
struct rhoUpdateParallel {
  using floatT = double;

  std::atomic<floatT> rho;
  floatT beta;
  UF *uf;

  rhoUpdateParallel(floatT betaa, UF *uff) :
    beta(betaa), uf(uff) {
    rho = std::numeric_limits<floatT>::max();}

  void run(nodeT *u, nodeT *v) {
    floatT myDist = nodeDistance(u, v);
    parlay::write_min(&rho, myDist, std::less<floatT>()); //check
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) return false; // filtering
    if (u->size()+v->size() <= beta) return false; // not in E_u, not considered for rho
    if (nodeDistance(u, v) >= rho) return false; // no subsequent finds can update rho
    return true;
  }

  bool start(nodeT *u) {
    if (u->size() > beta) {
      return true;
    } else {
      return false;// if node size < beta, so would children
    }
  }

  floatT getRho() { return rho;}

  bool wellSeparated(nodeT* u, nodeT* v, floatT s) {return geomWellSeparated(u, v, s);}
};

template<class nodeT, class UF>
struct wspGetParallel {
  using floatT = double;
  using bcpT = std::tuple<typename nodeT::objT*,
			  typename nodeT::objT*,
			  typename nodeT::objT::floatT>;
  using bufT = parBuf<bcpT>;

  floatT rhoLo;
  floatT rhoHi;
  floatT beta;
  nodeT *tree;
  UF *uf;
  bufT **out;

  wspGetParallel(floatT betaa, floatT rhoLoo, floatT rhoHii, nodeT *treee, UF *uff) :
    beta(betaa), rhoLo(rhoLoo), rhoHi(rhoHii), tree(treee), uf(uff) {
    size_t procs = parlay::num_workers();
    out = (bufT**) malloc(sizeof(bufT*)*procs);
    parlay::parallel_for(0, procs, [&](size_t p) {
			     out[p] = new bufT(tree->size()/procs);
			   });
  }

  ~wspGetParallel() {
    size_t procs = parlay::num_workers();
    parlay::parallel_for(0, procs, [&](size_t p) {
			     delete out[p];});
    free(out);
  }

  parlay::sequence<bcpT> collect() {
    int procs = parlay::num_workers();
      return parBufCollect<bcpT>(out, procs);
  }

  void run(nodeT *u, nodeT *v) {
    auto bcp = bichromaticClosestPair(u, v);
    if (u->size() + v->size() <= beta &&
	std::get<2>(bcp) >= rhoLo &&
        std::get<2>(bcp) < rhoHi) {
      auto tmp = out[parlay::worker_id()]->increment();
      get<0>(*tmp) = get<0>(bcp);
      get<1>(*tmp) = get<1>(bcp);
      get<2>(*tmp) = get<2>(bcp);
    }
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) {return false;}
    if (nodeDistance(u, v) >= rhoHi) return false; // too separated to be considered
    if (nodeFarDistance(u, v) < rhoLo) return false; // too close to be considered, bug!!
    return true;
  }

  bool start(nodeT *u) {
    if (u->diag() >= rhoLo) {
      return true;
    } else {
      return false;
    }
  }
  bool wellSeparated(nodeT* u, nodeT* v, floatT s) {return geomWellSeparated(u, v, s);}
};

template <class nodeT, class UF>
parlay::sequence<std::tuple<typename nodeT::objT*,
			    typename nodeT::objT*,
			    typename nodeT::objT::floatT>>
filterWspdParallel(double t_beta,
		   double t_rho_lo,
		   double& t_rho_hi,
		   nodeT *t_kdTree,
		   UF *t_mst) {
  using floatT = double;
  using objT = typename nodeT::objT;
  using bcpT = std::tuple<objT*, objT*, floatT>;

  auto myRho = rhoUpdateParallel<nodeT, UF>(t_beta, t_mst);

  pargeo::kdTree::computeWspdParallel<nodeT, rhoUpdateParallel<nodeT, UF>>(t_kdTree, &myRho);

  auto mySplitter = wspGetParallel<nodeT, UF>(t_beta, t_rho_lo, myRho.getRho(), t_kdTree, t_mst);

  pargeo::kdTree::computeWspdParallel<nodeT, wspGetParallel<nodeT, UF>>(t_kdTree, &mySplitter);

  t_rho_hi = myRho.getRho();
  return mySplitter.collect();
}

} // End namespace emstInternal
} // End namespace pargeo
