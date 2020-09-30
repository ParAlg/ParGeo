// This code is part of the project "Fast Parallel Algorithms for
// Euclidean Minimum SpanningTree and Hierarchical Spatial Clustering"
// Copyright (c) 2020 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
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

#ifndef UNREACHABLE_DECOMP
#define UNREACHABLE_DECOMP

#include <vector>
#include "kdTree.h"
#include "geometry.h"
#include "parBuf.h"
#include "wspd.h"

template<class nodeT, class objT>
inline void nodeCD(nodeT *nd, floatT *coreDist, floatT *cdMin, floatT *cdMax, nodeT* root, objT* P) {
  if(nd->isLeaf()){
    for (intT i = 0; i < nd->size(); ++i) {
      if (coreDist[nd->getItem(i)-P] > cdMax[nd-root]) {
        cdMax[nd-root] = coreDist[nd->getItem(i)-P];
      }
      if (coreDist[nd->getItem(i)-P] < cdMin[nd-root]) {
        cdMin[nd-root] = coreDist[nd->getItem(i)-P];
      }
    }
  } else {
    if (nd->size() > 2000) {
      cilk_spawn nodeCD(nd->L(), coreDist, cdMin, cdMax, root, P);
      nodeCD(nd->R(), coreDist, cdMin, cdMax, root, P);
      cilk_sync;
    } else {
      nodeCD(nd->L(), coreDist, cdMin, cdMax, root, P);
      nodeCD(nd->R(), coreDist, cdMin, cdMax, root, P);
    }
    cdMax[nd-root] = max(cdMax[nd->L()-root], cdMax[nd->R()-root]);
    cdMin[nd-root] = min(cdMin[nd->L()-root], cdMin[nd->R()-root]);
  }
}

//geo-separated || unreachable
template<int dim, class nodeT>
inline bool unreachable(nodeT *u, nodeT *v, floatT* cdMin, floatT* cdMax, nodeT *root) {
  floatT circleDiam_u = 0;
  floatT circleDiam_v = 0;
  floatT circleDistance = 0;
  for (int d = 0; d < dim; ++ d) {
    floatT uTmpDiff = u->getMax(d) - u->getMin(d);
    floatT vTmpDiff = v->getMax(d) - v->getMin(d);
    floatT uTmpAvg = (u->getMax(d) + u->getMin(d))/2;
    floatT vTmpAvg = (v->getMax(d) + v->getMin(d))/2;
    circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
    circleDiam_u += uTmpDiff * uTmpDiff;
    circleDiam_v += vTmpDiff * vTmpDiff;
  }
  circleDiam_u = sqrt(circleDiam_u);
  circleDiam_v = sqrt(circleDiam_v);

  floatT myRadius = max(circleDiam_u, circleDiam_v)/2;
  floatT myDiam = max(2*myRadius, cdMax[u-root]);
  myDiam = max(myDiam, cdMax[v-root]);

  circleDistance = sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;
  bool geoSep = circleDistance >= 2 * myRadius;
  circleDistance = max(circleDistance, cdMin[u-root]);
  circleDistance = max(circleDistance, cdMin[v-root]);

  if (circleDistance >= myDiam) {
    return true || geoSep;
  } else {
    return false || geoSep;
  }
}

template <int dim, class nodeT>
struct unreachableNormalSerial {
  typedef wsp<nodeT> pType;
  vector<pType> *out;
  floatT* cdMin;
  floatT* cdMax;
  nodeT* root;

  unreachableNormalSerial(nodeT* roott, vector<pType> *outt, floatT* cdMinn, floatT* cdMaxx) : root(roott), out(outt), cdMin(cdMinn), cdMax(cdMaxx) {}

  inline void run(nodeT *u, nodeT *v) {out->emplace_back(u, v);}
  inline bool moveon(nodeT *u, nodeT *v) {return true;}
  inline bool start(nodeT *u) {return true;}
  inline bool wellSeparated(nodeT *u, nodeT *v, int s) {return unreachable<dim, nodeT>(u, v, cdMin, cdMax, root);}
};

template <int dim, class nodeT>
struct unreachableNormalParallel {
  typedef wsp<nodeT> pType;
  typedef parBuf<pType> bufT;
  bufT **out;
  floatT* cdMin;
  floatT* cdMax;
  nodeT* root;

  unreachableNormalParallel(nodeT* roott, floatT* cdMinn, floatT* cdMaxx) : root(roott), cdMin(cdMinn), cdMax(cdMaxx) {
    int P = getWorkers();
    out = newA(bufT*, P);
    par_for(int p=0; p<P; ++p) {
      out[p] = new bufT(root->size()/P);
    }
  }

  vector<pType>* collect() {
    int P = getWorkers();
    auto tmp = parBufCollect<pType>(out, P);

    par_for(int p=0; p<P; ++p) delete out[p];
    free(out);
    return tmp;
  }

  inline void run(nodeT *u, nodeT *v) {
    auto tmp = out[getWorkerId()]->increment();
    tmp->u = u;
    tmp->v = v;
  }

  inline bool moveon(nodeT *u, nodeT *v) {return true;}
  inline bool start(nodeT *u) {return true;}
  inline bool wellSeparated(nodeT *u, nodeT *v, int s) {return unreachable<dim, nodeT>(u, v, cdMin, cdMax, root);}
};

#endif
