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

#ifndef KD_TREE_H
#define KD_TREE_H

#include <math.h>
#include <vector>
#include "pbbs/parallel.h"
#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "pbbs/parallel.h"
#include "geometry.h"

// *************************************************************
//   A generic parallel Euclidean kdTree
// *************************************************************
// objT needs to be templatized with int dim,
// and supports coordinate() that returns floatT[dim],
// and coordinate(i) that returns floatT* A[i]

template<int dim, class objT>
class kdNode {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef kdNode<dim, objT> nodeT;
  intT leafSize = 16;
  static const int boxInclude = 0;
  static const int boxOverlap = 1;
  static const int boxExclude = 2;
  int k;
  pointT pMin, pMax;
  objT **items;
  intT n;
  nodeT* left;
  nodeT* right;

  inline void boundingBoxSerial() {
    pMin = pointT(items[0]->coordinate());
    pMax = pointT(items[0]->coordinate());
    for(intT i=0; i<n; ++i) {
      pMin.minCoords(items[i]->coordinate());
      pMax.maxCoords(items[i]->coordinate());
    }}

  inline void boundingBoxParallel() {
    intT P = getWorkers()*8;
    intT blockSize = (n+P-1)/P;
    pointT localMin[P];
    pointT localMax[P];
    for (intT i=0; i<P; ++i) {
      localMin[i] = pointT(items[0]->coordinate());
      localMax[i] = pointT(items[0]->coordinate());}
    par_for(intT p=0; p<P; ++p) {
      intT s = p*blockSize;
      intT e = min((intT)(p+1)*blockSize,n);
      for (intT j=s; j<e; ++j) {
        localMin[p].minCoords(items[j]->coordinate());
        localMax[p].maxCoords(items[j]->coordinate());}
    }
    pMin = pointT(items[0]->coordinate());
    pMax = pointT(items[0]->coordinate());
    for(intT p=0; p<P; ++p) {
      pMin.minCoords(localMin[p].x);
      pMax.maxCoords(localMax[p].x);}
  }

  inline intT splitItemSerial(floatT xM) {
    if (n < 2) {
      cout << "error, kdTree splitting singleton, abort" << endl;abort();}
    intT lPt = 0;
    intT rPt = n-1;
    while (lPt < rPt) {
      if (items[lPt]->coordinate(k)>=xM) {
        while (items[rPt]->coordinate(k)>=xM && lPt < rPt) {
          rPt--;
        }
        if (lPt < rPt) {
          swap(items[lPt], items[rPt]);
          rPt--; }
        else { break;}
      }
      lPt++;
    }
    if (items[lPt]->coordinate(k) < xM) lPt++;
    return lPt;
  }

  inline intT splitItemParallel(floatT xM, objT **scratch, intT* flags) {
    if (n < 2) {
      cout << "error, kdTree splitting singleton, abort" << endl;abort();}
    par_for(intT i=0; i<n; ++i) {
      if (items[i]->coordinate(k)<xM) flags[i]=1;
      else flags[i] = 0;
    }
    intT leftSize = sequence::prefixSum(flags,0,n);
    par_for(intT i=0; i<n-1; ++i) {
      if (flags[i] != flags[i+1]) scratch[flags[i]] = items[i];
      if (i-flags[i] != i+1-flags[i+1]) scratch[leftSize+i-flags[i]] = items[i];
    }
    if (flags[n-1] != leftSize) scratch[flags[n-1]] = items[n-1];
    if (n-1-flags[n-1] != n-leftSize) scratch[leftSize+n-1-flags[n-1]] = items[n-1];
    par_for(intT i=0; i<n; ++i) {
      items[i] = scratch[i];
    }
    return leftSize;
  }

  inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2) {
    bool exclude = false;
    bool include = true;//1 include 2
    for(int i=0; i<dim; ++i) {
      if (pMax1[i]<pMin2[i] || pMin1[i]>pMax2[i]) exclude = true;
      if (pMax1[i]<pMax2[i] || pMin1[i]>pMin2[i]) include = false;
    }
    if (exclude) return boxExclude;
    else if (include) return boxInclude;
    else return boxOverlap;
  }

  inline bool itemInBox(pointT pMin1, pointT pMax1, objT* item) {
    for(int i=0; i<dim; ++i) {
      if (pMax1[i]<item->coordinate(i) || pMin1[i]>item->coordinate(i)) return false;
    }
    return true;}

  void constructSerial(nodeT *space, bool noCoarsen) {
    boundingBoxSerial();
    if (isLeaf()) {
      left = NULL; right = NULL;
    } else {
      if (!space[0].isEmpty() || !space[1].isEmpty()) {
        cout << "error, kdNode overwrite, abort" << endl;abort();}
      floatT xM = -1;
      for (int kk=0; kk<dim; ++kk) {
        if (pMax[kk]-pMin[kk]>xM) {
          xM = pMax[kk]-pMin[kk];
          k = kk;}}
      xM = (pMax[k]+pMin[k])/2;
      intT median = splitItemSerial(xM);
      if (median == 0 || median == n) {median = n/2;}
      space[0] = nodeT(items, median, space+1, noCoarsen);
      space[2*median-1] = nodeT(items+median, n-median, space+2*median, noCoarsen);
      left = space;
      right = space+2*median-1;
    }
  }

  //cilk_spawn requires function
  void buildNode(nodeT *space, objT** itemss, intT nn, nodeT *spacee, objT** scratchh, intT* flagss, bool noCoarsen) {
    space[0] = nodeT(itemss, nn, spacee, scratchh, flagss, noCoarsen);}
  void constructParallel(nodeT *space, objT** scratch, intT* flags, bool noCoarsen) {
    boundingBoxParallel();
    if (isLeaf()) {
      left = NULL; right = NULL;
    } else {
      if (!space[0].isEmpty() || !space[1].isEmpty()) {
        cout << "error, kdNode overwrite, abort" << endl;abort();}
      floatT xM = -1;
      for (int kk=0; kk<dim; ++kk) {
        if (pMax[kk]-pMin[kk]>xM) {
          xM = pMax[kk]-pMin[kk];
          k = kk;}}
      xM = (pMax[k]+pMin[k])/2;
      intT median = splitItemParallel(xM, scratch, flags);
      if (median == 0 || median == n) {median = n/2;}
      cilk_spawn buildNode(&space[0], items, median, space+1, scratch, flags, noCoarsen);
      buildNode(&space[2*median-1], items+median, n-median, space+2*median, scratch+median, flags+median, noCoarsen);
      cilk_sync;
      left = space;
      right = space+2*median-1;
    }
  }

  public:
  nodeT* L() {return left;}
  nodeT* R() {return right;}
  inline intT size() {return n;}
  inline objT* operator[](intT i) {return items[i];}

  kdNode(objT** itemss, intT nn, nodeT *space, objT** scratch, intT* flags, bool noCoarsen=false): items(itemss), n(nn) {
    if (noCoarsen) leafSize = 1;
    if (n>2000) constructParallel(space, scratch, flags, noCoarsen);
    else constructSerial(space, noCoarsen);}
  kdNode(objT** itemss, intT nn, nodeT *space, bool noCoarsen=false): items(itemss), n(nn) {
    if (noCoarsen) leafSize = 1;
    constructSerial(space, noCoarsen);}

  void setEmpty() {n=-1;}
  bool isEmpty() {return n<0;}
  bool isLeaf() {return n<=leafSize;}

  floatT nodeDiag() {//todo change name
    floatT result = 0;
    for (int d = 0; d < dim; ++ d) {
      floatT tmp = pMax[d] - pMin[d];
      result += tmp * tmp;
    }
    return sqrt(result);
  }

  // return maximum span of node bounding box among all dimensions
  inline floatT lMax() {
    floatT myMax = 0;
    for (int d=0; d<dim; ++d) {
      floatT thisMax = pMax[d] - pMin[d];
      if (thisMax > myMax) {
	myMax = thisMax;}
    }
    return myMax;
  }

  //return the bb distance with and n2
  inline floatT nodeDistance(nodeT* n2) {
    for (int d = 0; d < dim; ++ d) {
      if (pMin[d] > n2->pMax[d] || n2->pMin[d] > pMax[d]) {
        // disjoint at dim d, and intersect on dim < d
        floatT rsqr = 0;
        for (int dd = d; dd < dim; ++ dd) {
          floatT tmp = max(pMin[dd]-n2->pMax[dd], n2->pMin[dd]-pMax[dd]);
          tmp = max(tmp, (floatT)0);
          rsqr += tmp*tmp;
        }
        return sqrt(rsqr);
      }
    }
    return 0; // intersect
  }

  // inline floatT objDistance(objT* p2) {
  //   for (int d = 0; d < dim; ++ d) {
  //     if (pMin[d] > p2->coordinate(d) || p2->coordinate(d) > pMax[d]) {
  //       floatT rsqr = 0;
  //       for (int dd = d; dd < dim; ++ dd) {
  //         floatT tmp = max(pMin[dd]-p2->coordinate(dd), p2->coordinate(dd)-pMax[dd]);
  //         tmp = max(tmp, (floatT)0);
  //         rsqr += tmp*tmp;
  //       }
  //       return sqrt(rsqr);
  //     }
  //   }
  //   return 0; // p2 inside
  // }

  //return the far bb distance between n1 and n2
  inline floatT nodeFarDistance(nodeT* n2) {
    floatT result = 0;
    for (int d = 0; d < dim; ++ d) {
      floatT tmp = max(pMax[d],n2->pMax[d]) - min(pMax[d],n2->pMax[d]);
      result += tmp *tmp;
    }
    return sqrt(result);
  }

  //whether well separated with v
  inline bool wellSeparated(nodeT *v) {
    static const int s = 2; // separation constant
    floatT circleDiam_u = 0;
    floatT circleDiam_v = 0;
    floatT circleDistance = 0;
    for (int d = 0; d < dim; ++ d) {
      floatT uTmpDiff = pMax[d] - pMin[d];
      floatT vTmpDiff = v->pMax[d] - v->pMin[d];
      floatT uTmpAvg = (pMax[d] + pMin[d])/2;
      floatT vTmpAvg = (v->pMax[d] + v->pMin[d])/2;
      circleDistance += (uTmpAvg - vTmpAvg) * (uTmpAvg - vTmpAvg);
      circleDiam_u += uTmpDiff * uTmpDiff;
      circleDiam_v += vTmpDiff * vTmpDiff;
    }
    circleDiam_u = sqrt(circleDiam_u);
    circleDiam_v = sqrt(circleDiam_v);
    floatT myRadius = max(circleDiam_u, circleDiam_v)/2;
    circleDistance = sqrt(circleDistance) - circleDiam_u/2 - circleDiam_v/2;
    return circleDistance >= (s * myRadius);
  }

  //Deprecate
  void rangeNeighbor(pointT pMin1, pointT pMax1, floatT r, vector<objT*>* accum) {
    int relation = boxCompare(pMin1, pMax1, pMin, pMax);
    if (relation == boxInclude) {
      for(intT i=0; i<n; ++i) accum->push_back(items[i]);
    } else if (relation == boxOverlap) {
      if (isLeaf()) {
        for(intT i=0; i<n; ++i) {
          if (itemInBox(pMin1, pMax1, items[i])) accum->push_back(items[i]);}
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, accum);
        right->rangeNeighbor(pMin1, pMax1, r, accum);}
    }
  }

  //Deprecate
  template<class func>
  void rangeNeighbor(pointT pMin1, pointT pMax1, floatT r, func* f) {
    if (f->isComplete()) return;
    int relation = boxCompare(pMin1, pMax1, pMin, pMax);
    if (relation == boxInclude) {
      for(intT i=0; i<n; ++i) {
        if (f->checkComplete(items[i])) break;
      }
    } else if (relation == boxOverlap) {
      if (isLeaf()) {
        for(intT i=0; i<n; ++i) {
          if (itemInBox(pMin1, pMax1, items[i])) {
            if (f->checkComplete(items[i])) break;}
        }
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, f);
        right->rangeNeighbor(pMin1, pMax1, r, f);}
    }
  }

  template<class func, class func2>
  void rangeNeighbor(pointT pMin1, pointT pMax1, floatT r, func term, func2 doTerm) {
    if (term()) return;
    int relation = boxCompare(pMin1, pMax1, pMin, pMax);
    if (relation == boxInclude) {
      for(intT i=0; i<n; ++i) {
        if (doTerm(items[i])) break;
      }
    } else if (relation == boxOverlap) {
      if (isLeaf()) {
        for(intT i=0; i<n; ++i) {
          if (itemInBox(pMin1, pMax1, items[i])) {
            if (doTerm(items[i])) break;}
        }
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, term, doTerm);
        right->rangeNeighbor(pMin1, pMax1, r, term, doTerm);}
    }
  }

  struct bcp {
    objT* u;
    objT* v;
    floatT dist;
    bcp(objT* uu, objT* vv, floatT distt): u(uu), v(vv), dist(distt) {}
    bcp(): u(NULL), v(NULL), dist(floatMax()) {}
    void update(objT* uu, objT* vv) {
      auto distt = uu->dist(*vv);
      if (distt < dist) {
        u = uu; v = vv; dist = distt;}
    }
  };

  inline void compBcpH(nodeT* n2, bcp* r) {
    if (nodeDistance(n2) > r->dist) return;

    if (isLeaf() && n2->isLeaf()) {//basecase
      for (intT i=0; i<size(); ++i) {
        for (intT j=0; j<n2->size(); ++j) {
          r->update(items[i], n2->items[j]);}
      }
    } else {//recursive, todo consider call order, might help
      if (isLeaf()) {
        compBcpH(n2->left, r);
        compBcpH(n2->right, r);
      } else if (n2->isLeaf()) {
        n2->compBcpH(left, r);
        n2->compBcpH(right, r);
      } else {
        left->compBcpH(n2->left, r);
        left->compBcpH(n2->right, r);
        right->compBcpH(n2->left, r);
        right->compBcpH(n2->right, r);}
    }
  }

  inline bcp compBcp(nodeT* n2) {
    auto r = bcp();
    compBcpH(n2, &r);
    // for (intT i=0; i<size(); ++i) {
    //   for (intT j=0; j<n2->size(); ++j) {
    //     r.update(items[i], n2->items[j]);}
    // }
    return r;
  }

};

template<int dim, class objT>
class kdTree {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef kdNode<dim, objT> nodeT;
  objT **items;
  nodeT *root;

  public:
  nodeT* rootNode() {return root;}

  kdTree(objT* P, intT n, bool parallel=true, bool noCoarsen=false) {
    items = newA(objT*, n);
    par_for(intT i=0; i<n; ++i) {
      items[i]=&P[i];
    }
    root = newA(nodeT, 2*n-1);
    par_for(intT i=0; i<2*n-1; ++i) {
      root[i].setEmpty();
    }
    if (parallel) {
      objT** scratch = newA(objT*, n);
      intT* flags = newA(intT, n);
      root[0] = nodeT(items, n, root+1, scratch, flags, noCoarsen);
      free(scratch);
      free(flags);
    } else {
      root[0] = nodeT(items, n, root+1, noCoarsen);}
  }
  ~kdTree() {
    free(items);
    free(root);}

  vector<objT*>* rangeNeighbor(objT* query, floatT r) {
    vector<objT*>* accum = new vector<objT*>();
    pointT pMin1 = pointT();
    pointT pMax1 = pointT();
    floatT* center = query->coordinate();
    for (int i=0; i<dim; ++i) {
      pMin1.updateX(i, center[i]-r);
      pMax1.updateX(i, center[i]+r);}
    root->rangeNeighbor(pMin1, pMax1, r, accum);
    return accum;
  }

  //Deprecate
  template<class func>
  void rangeNeighbor(objT* query, floatT r, func* f) {
    pointT pMin1 = pointT();
    pointT pMax1 = pointT();
    floatT* center = query->coordinate();
    for (int i=0; i<dim; ++i) {
      pMin1.updateX(i, center[i]-r);
      pMax1.updateX(i, center[i]+r);}
    root->rangeNeighbor(pMin1, pMax1, r, f);
  }

  template<class func, class func2>
  void rangeNeighbor(objT* query, floatT r, func term, func2 doTerm) {
    pointT pMin1 = pointT();
    pointT pMax1 = pointT();
    floatT* center = query->coordinate();
    for (int i=0; i<dim; ++i) {
      pMin1.updateX(i, center[i]-r);
      pMax1.updateX(i, center[i]+r);}
    root->rangeNeighbor(pMin1, pMax1, r, term, doTerm);
  }
  // objT* nearestNeighbor(objT*) {};
  // vector<objT*> kNearestNeighbor(objT*, intT k) {};
};

#endif
