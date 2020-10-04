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

#ifndef BCCP_STAR_H
#define BCCP_STAR_H

//computes the bccp* in terms of mutual reachability

template<class nodeT, class bcpT, class objT>
inline void compBcpStarH(nodeT* n1, nodeT* n2, bcpT* r, floatT* coreDist, objT* P) {
  if (n1->nodeDistance(n2) > r->dist) return;

  if (n1->isLeaf() && n2->isLeaf()) {//basecase
    for (intT i=0; i<n1->size(); ++i) {
      for (intT j=0; j<n2->size(); ++j) {
        floatT dist = max(n1->getItem(i)->dist(*n2->getItem(j)),
                          coreDist[n1->getItem(i)-P]);
        dist = max(dist, coreDist[n2->getItem(j)-P]);
        r->update(n1->getItem(i), n2->getItem(j), dist);
      }
    }
  } else {//recursive, todo consider call order, might help
    if (n1->isLeaf()) {
      if (n1->nodeDistance(n2->L()) < n1->nodeDistance(n2->R())) {
        compBcpStarH(n1, n2->L(), r, coreDist, P);
        compBcpStarH(n1, n2->R(), r, coreDist, P);
      } else {
        compBcpStarH(n1, n2->R(), r, coreDist, P);
        compBcpStarH(n1, n2->L(), r, coreDist, P);
      }
    } else if (n2->isLeaf()) {
      if (n2->nodeDistance(n1->L()) < n2->nodeDistance(n1->R())) {
        compBcpStarH(n2, n1->L(), r, coreDist, P);
        compBcpStarH(n2, n1->R(), r, coreDist, P);
      } else {
        compBcpStarH(n2, n1->R(), r, coreDist, P);
        compBcpStarH(n2, n1->L(), r, coreDist, P);
      }
    } else {
      pair<nodeT*, nodeT*> ordering[4];
      ordering[0] = make_pair(n2->L(), n1->L());
      ordering[1] = make_pair(n2->R(), n1->L());
      ordering[2] = make_pair(n2->L(), n1->R());
      ordering[3] = make_pair(n2->R(), n1->R());
      auto bbd = [&](pair<nodeT*,nodeT*> p1, pair<nodeT*,nodeT*> p2) {
                                                                      return p1.first->nodeDistance(p1.second) < p2.first->nodeDistance(p2.second);};
      quickSortSerial(ordering, 4, bbd);
      for (intT o=0; o<4; ++o) {
        compBcpStarH(ordering[o].first, ordering[o].second, r, coreDist, P);}
    }
  }
}

template<class nodeT, class bcpT, class objT>
inline bcpT compBcpStarBF(nodeT* n1, nodeT* n2, floatT* coreDist, objT* P) {
  auto r = bcpT();
  for (intT i=0; i<n1->size(); ++i) {
    for (intT j=0; j<n2->size(); ++j) {
      floatT dist = max(n1->getItem(i)->dist(*n2->getItem(j)),
                        coreDist[n1->getItem(i)-P]);
      dist = max(dist, coreDist[n2->getItem(j)-P]);
      r.update(n1->getItem(i), n2->getItem(j), dist);
    }
  }
  return r;
}

template<class nodeT, class bcpT, class objT>
inline bcpT compBcpStar(nodeT* n1, nodeT* n2, floatT* coreDist, objT* P) {
  if (n1->size() + n2->size() <= 32) {
    return compBcpStarBF<nodeT, bcpT, objT>(n1, n2, coreDist, P);}

  auto r = bcpT();
  compBcpStarH(n1, n2, &r, coreDist, P);
  return r;
}

#endif
