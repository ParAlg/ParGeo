// Copyright (c) 2020 Yiqiu Wang
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

#ifndef WSPD_FILTER
#define WSPD_FILTER

#include <vector>
#include "parBuf.h"
#include "kdTree.h"
#include "kdNode.h"
#include "pbbs/unionFind.h"
#include "pbbs/utils.h"

template<class nodeT>
struct rhoUpdateSerial {
  floatT rho;
  floatT beta;
  edgeUnionFind *uf;

  rhoUpdateSerial(floatT betaa, edgeUnionFind *uff) :
    beta(betaa), uf(uff) {rho = floatMax();}

  void run(nodeT *u, nodeT *v) {
    floatT myDist = u->nodeDistance(v);
    if (myDist < rho) rho = myDist;
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) return false; // filtering
    if (u->size()+v->size() <= beta) return false; // not in E_u, not considered for rho
    if (u->nodeDistance(v) >= rho) return false; // no subsequent finds can update rho
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
  bool wellSeparated(nodeT* u, nodeT* v, int s) {return u->wellSeparated(v,s);}
};

template<class treeT, class nodeT>
struct wspGetSerial {
  typedef struct nodeT::bcp bcpT;

  floatT rhoLo;
  floatT rhoHi;
  floatT beta;
  vector<bcpT> *out;
  treeT *tree;
  edgeUnionFind *uf;

  wspGetSerial(floatT betaa, floatT rhoLoo, floatT rhoHii, vector<bcpT> *outt, treeT *treee, edgeUnionFind *uff) :
  beta(betaa), rhoLo(rhoLoo), rhoHi(rhoHii), out(outt), tree(treee), uf(uff) {}

  void run(nodeT *u, nodeT *v) {
    auto bcp = u->compBcp(v);
    if (u->size() + v->size() <= beta &&
        bcp.dist >= rhoLo &&
        bcp.dist < rhoHi) {
      out->push_back(bcp);
    }
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) {return false;}
    if (u->nodeDistance(v) >= rhoHi) return false; // too separated to be considered
    if (u->nodeFarDistance(v) < rhoLo) return false; // too close to be considered, bug!!
    return true;
  }

  bool start(nodeT *u) {
    if (u->nodeDiag() >= rhoLo) {
      return true;
    } else {
      return false;
    }
  }
  bool wellSeparated(nodeT* u, nodeT* v, int s) {return u->wellSeparated(v,s);}
};

template <class treeT, class nodeT>
inline floatT filterWspdSerial(floatT t_beta, floatT t_rho_lo, vector<typename nodeT::bcp> *t_wsp, treeT *t_kdTree, edgeUnionFind *t_mst) {
  auto myRho = rhoUpdateSerial<nodeT>(t_beta, t_mst);
  wspdSerial<nodeT, rhoUpdateSerial<nodeT>>(t_kdTree->rootNode(), &myRho);
  auto mySplitter = wspGetSerial<treeT, nodeT>(t_beta, t_rho_lo, myRho.getRho(), t_wsp, t_kdTree, t_mst);
  wspdSerial<nodeT, wspGetSerial<treeT, nodeT>>(t_kdTree->rootNode(), &mySplitter);
  return myRho.getRho();
}


template<class nodeT>
struct rhoUpdateParallel {
  floatT rho;
  floatT beta;
  edgeUnionFind *uf;

  rhoUpdateParallel(floatT betaa, edgeUnionFind *uff) :
    beta(betaa), uf(uff) {rho = floatMax();}

  void run(nodeT *u, nodeT *v) {
    floatT myDist = u->nodeDistance(v);
    utils::writeMin(&rho, myDist);
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) return false; // filtering
    if (u->size()+v->size() <= beta) return false; // not in E_u, not considered for rho
    if (u->nodeDistance(v) >= rho) return false; // no subsequent finds can update rho
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
  bool wellSeparated(nodeT* u, nodeT* v, int s) {return u->wellSeparated(v,s);}
};

template<class treeT, class nodeT>
struct wspGetParallel {
  typedef struct nodeT::bcp bcpT;
  typedef parBuf<bcpT> bufT;

  floatT rhoLo;
  floatT rhoHi;
  floatT beta;
  bufT **out;
  treeT *tree;
  edgeUnionFind *uf;

  wspGetParallel(floatT betaa, floatT rhoLoo, floatT rhoHii, treeT *treee, edgeUnionFind *uff) :
    beta(betaa), rhoLo(rhoLoo), rhoHi(rhoHii), tree(treee), uf(uff) {
    int P = getWorkers();
    out = newA(bufT*, P);
    intT n = tree->rootNode()->size();
    par_for(int p=0; p<P; ++p) {
      out[p] = new bufT(n/P);
    }
  }

  vector<bcpT>* collect() {
    int P = getWorkers();
    auto tmp = parBufCollect<bcpT>(out, P);

    par_for(int p=0; p<P; ++p) delete out[p];
    free(out);
    return tmp;
  }

  void run(nodeT *u, nodeT *v) {
    auto bcp = u->compBcp(v);
    if (u->size() + v->size() <= beta &&
        bcp.dist >= rhoLo &&
        bcp.dist < rhoHi) {
      auto tmp = out[getWorkerId()]->increment();
      tmp->u = bcp.u; tmp->v = bcp.v; tmp->dist = bcp.dist;
    }
  }

  bool moveon(nodeT *u, nodeT *v) {
    if (u->hasId() && u->getId() == v->getId()) {return false;}
    if (u->nodeDistance(v) >= rhoHi) return false; // too separated to be considered
    if (u->nodeFarDistance(v) < rhoLo) return false; // too close to be considered, bug!!
    return true;
  }

  bool start(nodeT *u) {
    if (u->nodeDiag() >= rhoLo) {
      return true;
    } else {
      return false;
    }
  }
  bool wellSeparated(nodeT* u, nodeT* v, int s) {return u->wellSeparated(v,s);}
};

template <class treeT, class nodeT>
inline pair<floatT, vector<typename nodeT::bcp>*> filterWspdParallel(floatT t_beta, floatT t_rho_lo, treeT *t_kdTree, edgeUnionFind *t_mst) {
  auto myRho = rhoUpdateParallel<nodeT>(t_beta, t_mst);
  wspdParallel<nodeT, rhoUpdateParallel<nodeT>>(t_kdTree->rootNode(), &myRho);
  auto mySplitter = wspGetParallel<treeT, nodeT>(t_beta, t_rho_lo, myRho.getRho(), t_kdTree, t_mst);
  wspdParallel<nodeT, wspGetParallel<treeT, nodeT>>(t_kdTree->rootNode(), &mySplitter);
  return make_pair(myRho.getRho(), mySplitter.collect());
}

#endif
