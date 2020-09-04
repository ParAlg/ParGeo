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

#ifndef DYNAMIC_KD_TREE_H
#define DYNAMIC_KD_TREE_H

#include <math.h>
#include <vector>
#include "geometry.h"
#include "pbbs/utils.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/parallel.h"

// *************************************************************
//   A generic parallel dynamic Euclidean kdTree
// *************************************************************
// objT needs to be templatized with int dim,
// and supports coordinate() that returns floatT[dim],
// and coordinate(i) that returns floatT* A[i]

template<int dim> struct geoPointHash;

template<int dim, class objT>
class dynamicKdNode {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef dynamicKdNode<dim, objT> nodeT;
  struct splitter {
    int k; double x;
    splitter(): k(-1), x(-1) {}
    splitter(int kk, double xx): k(kk), x(xx) {}
    inline bool leftOf(objT* item) {return item->coordinate(k) < x;}
    inline bool rightOf(objT* item) {return !leftOf(item);}
  };
  typedef splitter splitterT;
  static const intT leafSize = 32;
  static const int boxInclude = 0;
  static const int boxOverlap = 1;
  static const int boxExclude = 2;
  int k;
  pointT pMin, pMax;
  objT **items;
  intT numItems;
  intT numErased;
  intT capacity;
  nodeT* left;
  nodeT* right;
  splitterT split;
  intT counter;

  inline void boundingBoxExpandSerial(objT** newItems, intT nn) {
    for(intT i=0; i<nn; ++i) {
      if (newItems[i]) {
        pMin.minCoords(newItems[i]->coordinate());
        pMax.maxCoords(newItems[i]->coordinate());}
    }}

  inline void boundingBoxExpandParallel(objT** newItems, intT nn) {
    intT P = getWorkers()*8;
    intT blockSize = (nn+P-1)/P;
    pointT localMin[P];
    pointT localMax[P];
    intT first = 0;
    while(!newItems[first] && first < nn) first++;
    if(first >= nn) return;
    for (intT i=0; i<P; ++i) {
      localMin[i] = pointT(newItems[first]->coordinate());
      localMax[i] = pointT(newItems[first]->coordinate());}
    par_for(intT p=0; p<P; ++p) {
      intT s = p*blockSize;
      intT e = min((intT)(p+1)*blockSize,nn);
      for (intT j=s; j<e; ++j) {
        if (newItems[j]) {
          localMin[p].minCoords(newItems[j]->coordinate());
          localMax[p].maxCoords(newItems[j]->coordinate());}
      }
    }
    for(intT p=0; p<P; ++p) {
      pMin.minCoords(localMin[p].x);
      pMax.maxCoords(localMax[p].x);}
  }

  inline intT splitItemSerial(objT** itemss, intT nn, splitterT split) {
    if (nn < 2) {
      if (!itemss[0]) {
        return 0;//arbitrary
      } else if (itemss[0]->coordinate(split.k) < split.x) {
        return 1;
      } else return 0;
    }
    intT lPt = 0;
    intT rPt = nn-1;
    while (lPt < rPt) {
      if (itemss[lPt] && itemss[lPt]->coordinate(split.k)>=split.x) {
        while (itemss[rPt] && itemss[rPt]->coordinate(split.k)>=split.x && lPt < rPt) {
          rPt--;
        }
        if (lPt < rPt) {
          swap(itemss[lPt], itemss[rPt]);
          rPt--; }
        else { break;}
      }
      lPt++;
    }
    if (itemss[lPt] && itemss[lPt]->coordinate(split.k) < split.x) lPt++;
    return lPt;
  }

  inline intT splitItemParallel(objT** itemss, intT nn, splitterT split, objT **scratch=NULL, intT* flags=NULL) {
    if (nn < 2) {
      if (!itemss[0]) {
        return 0;//arbitrary
      } if (itemss[0]->coordinate(split.k) < split.x) {
        return 1;
      } else return 0;
    }

    if (nn < 2000) return splitItemSerial(itemss, nn, split);

    bool freeScratch = false;
    if(!scratch) {
      scratch = newA(objT*, nn);
      freeScratch = true;}
    bool freeFlag = false;
    if(!scratch) {
      flags = newA(intT, nn);
      freeFlag = true;}

    par_for(intT i=0; i<nn; ++i) {
      if (itemss[i] && itemss[i]->coordinate(split.k)<split.x) flags[i]=1;
      else flags[i] = 0;
    }
    intT leftSize = sequence::prefixSum(flags,0,nn);
    par_for(intT i=0; i<nn-1; ++i) {
      if (flags[i] != flags[i+1]) scratch[flags[i]] = itemss[i];
      if (i-flags[i] != i+1-flags[i+1]) scratch[leftSize+i-flags[i]] = itemss[i];
    }
    if (flags[nn-1] != leftSize) scratch[flags[nn-1]] = itemss[nn-1];
    if (nn-1-flags[nn-1] != nn-leftSize) scratch[leftSize+nn-1-flags[nn-1]] = itemss[nn-1];
    par_for(intT i=0; i<nn; ++i) {
      itemss[i] = scratch[i];
    }

    if(freeScratch) free(scratch);
    if(freeFlag) free(flags);

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

  void constructSerial(objT** itemss, intT nn) {
    pMin = pointT(itemss[0]->coordinate());
    pMax = pointT(itemss[0]->coordinate());
    boundingBoxExpandSerial(itemss, nn);
    if (nn <= leafSize) {
      resizeSerial(nn*2);
      for(intT i=0; i<nn; ++i) {items[i]=itemss[i];}
      numItems = nn;
      left = NULL; right = NULL;
    } else {
      floatT xM = -1;
      for (int kk=0; kk<dim; ++kk) {
        if (pMax[kk]-pMin[kk]>xM) {
          xM = pMax[kk]-pMin[kk];
          k = kk;}}
      xM = (pMax[k]+pMin[k])/2;
      split = splitterT(k, xM);
      intT median = splitItemSerial(itemss, nn, split);
      if(median == 0 || median == nn) {
        leaf = true;
        resizeSerial(nn*2);
        numItems = nn;
        for(intT i=0; i<nn; ++i) {items[i]=itemss[i];}
        left = NULL; right = NULL;
      } else {
        leaf = false;
        numItems = nn;
        left = new nodeT(itemss, median);
        right = new nodeT(itemss+median, nn-median);}
    }
  }

  //cilk_spawn requires function
  void buildLeftNode(objT** itemss, intT nn, objT** scratchh, intT* flagss) {
    left = new nodeT(itemss, nn, scratchh, flagss);}
  void buildRightNode(objT** itemss, intT nn, objT** scratchh, intT* flagss) {
    right = new nodeT(itemss, nn, scratchh, flagss);}
  void constructParallel(objT** itemss, intT nn, objT** scratch, intT* flags) {
    pMin = pointT(itemss[0]->coordinate());
    pMax = pointT(itemss[0]->coordinate());
    boundingBoxExpandParallel(itemss, nn);
    if (nn <= leafSize) {
      resizeSerial(nn*2);
      for(intT i=0; i<nn; ++i) {items[i]=itemss[i];}
      numItems = nn;
      left = NULL; right = NULL;
    } else {
      floatT xM = -1;
      for (int kk=0; kk<dim; ++kk) {
        if (pMax[kk]-pMin[kk]>xM) {
          xM = pMax[kk]-pMin[kk];
          k = kk;}}
      xM = (pMax[k]+pMin[k])/2;
      split = splitterT(k, xM);
      intT median = splitItemParallel(itemss, nn, split, scratch, flags);
      if(median == 0 || median == nn) {
        leaf = true;
        resizeSerial(nn*2);
        numItems = nn;
        for(intT i=0; i<nn; ++i) {items[i]=itemss[i];}
        left = NULL; right = NULL;
      } else {
        leaf = false;
        numItems = nn;
        cilk_spawn buildLeftNode(itemss, median, scratch, flags);
        buildRightNode(itemss+median, nn-median, scratch+median, flags+median);
        cilk_sync;
      }
    }
  }

  void compactSerial() {
    intT ii=0;
    for(intT i=0; i < numItems; ++i) {
      if(items[i]) items[ii++]=items[i];}
    numItems = ii;
    numErased = 0;
    for(intT i=numItems; i<capacity; ++i) items[i]=NULL;
  }

  void compactParallel(objT** scratch=NULL, intT* flags=NULL) {
    if (numItems < 2000) return compactSerial();

    bool freeScratch=false;
    if(!scratch) {
      scratch = newA(objT*, numItems);
      freeScratch = true;}
    bool freeFlag=false;
    if(!flags) {
      flags = newA(intT, numItems);
      freeFlag = true;}

    par_for(intT i=0; i<numItems; ++i) {scratch[i] = items[i];}
    par_for(intT i=0; i<numItems; ++i) {
      if (items[i]) flags[i] = 1;
      else flags[i] = 0;
    }
    intT valid = sequence::prefixSum(flags,0,numItems);
    par_for(intT i=0; i<numItems-1; ++i) {
      if (flags[i]!=flags[i+1]) {
        items[flags[i]] = scratch[i];}
    }
    intT i=numItems-1;
    if (flags[i]!=valid) items[flags[i]] = scratch[i];
    par_for(intT i=valid; i<capacity; ++i) {
      items[i] = NULL;
    }
    numItems = valid;
    numErased = 0;

    if(freeScratch) free(scratch);
    if(freeFlag) free(flags);
  }

  //free buffer, but still keep track of descendants size and bbox
  void makeInternal() {
    numItems = numItems-numErased;
    numErased = 0;
    capacity = 0;
    leaf = false;
    free(items);}

  void resizeSerial(intT want) {
    auto oldCapacity = capacity;
    capacity = max(want, numItems);
    if (capacity <= 0) capacity = 10;
    auto newItems = newA(objT*, capacity);
    for(intT i=0; i < numItems; ++i) newItems[i]=items[i];
    for(intT i=numItems; i < capacity; ++i) newItems[i]=NULL;
    swap(newItems, items);
    if (oldCapacity>0) free(newItems);
  }

  void resizeParallel(intT want) {
    auto oldCapacity = capacity;
    capacity = max(want, numItems);
    if (capacity <= 0) capacity = 10;
    auto newItems = newA(objT*, capacity);
    par_for(intT i=0; i<numItems; ++i) {newItems[i]=items[i];}
    par_for(intT i=numItems; i<capacity; ++i) {newItems[i]=NULL;}
    swap(newItems, items);
    if (oldCapacity>0) free(newItems);
  }

  //collect from all descendants
  void collectItems(objT** buffer) {
    if(!isLeaf()) {
      if(actualSize() > 2000) {
        cilk_spawn left->collectItems(buffer);
        right->collectItems(buffer+left->actualSize());
        cilk_sync;
      } else {
        left->collectItems(buffer);
        right->collectItems(buffer+left->actualSize());
      }
    } else {
      intT ii=0;
      for(intT i=0; i<size(); ++i) {
        if(items[i]){
          buffer[ii++]=items[i];
        }
      }
    }
  }

  void freeDescendant() {
    if(!isLeaf()) {
      if(left) {
        left->freeDescendant();
        delete left;
      }
      if(right) {
        right->freeDescendant();
        delete right;
      }
    }
  }

  void rebuild() {
    intT nn = actualSize();
    if(nn <= 0) return;
    objT** itemss = newA(objT*, nn*2);
    objT** scratch = itemss+nn;
    intT* flags = newA(intT, nn);
    collectItems(itemss);

    //reinitialize
    freeDescendant();
    capacity = 0;
    numItems = 0;
    numErased = 0;
    counter = 0;
    items = NULL;
    leaf = true;

    if(actualSize() > 2000) constructParallel(itemss, actualSize(), scratch, flags);
    else constructSerial(itemss, nn);
    free(itemss);
    free(flags);
  }

  public:
  dynamicKdNode(objT** itemss, intT nn, objT** scratch, intT* flags): capacity(0), numItems(0), numErased(0), leaf(true), counter(0), items(NULL) {
    if (nn>2000) constructParallel(itemss, nn, scratch, flags);
    else constructSerial(itemss, nn);
  }
  dynamicKdNode(objT** itemss, intT nn): capacity(0), numItems(0), numErased(0), leaf(true), counter(0), items(NULL) {
    constructSerial(itemss, nn);}
  ~dynamicKdNode() {
    if(items && isLeaf()) free(items);
  }

  bool leaf;
  inline void setEmpty() {capacity=-1;}
  inline bool isEmpty() {return capacity<0;}
  inline bool isLeaf() {return leaf;}
  inline intT actualSize() {return numItems-numErased;}
  inline intT size() {return numItems;}
  inline nodeT* leftChild() {return left;}
  inline nodeT* rightChild() {return right;}

  void rebuildRecurse() {
    if (floatT(counter)/floatT(actualSize()) >= 0.75) {
      rebuild();
    } else if(!isLeaf()) {
      if(actualSize()<=2000) {
        left->rebuildRecurse();
        right->rebuildRecurse();
      } else {
        cilk_spawn left->rebuildRecurse();
        right->rebuildRecurse();
        cilk_sync;
      }
    }
  }

  void erase(objT** itemss, intT nn) {
    if(nn <= 0) return;
    counter += nn;
    if(isLeaf()) {//small
      for(intT j=0; j<nn; ++j) {
        for(intT i=0; i<numItems; ++i) {
          if(items[i] && items[i] == itemss[j]) {
            items[i]=NULL;
            numErased++;}
        }}
    } else {
      intT median = splitItemSerial(itemss, nn, split);
      if(nn < 2000) {
        left->erase(itemss, median);
        right->erase(itemss+median, nn-median);
      } else {
        cilk_spawn left->erase(itemss, median);
        right->erase(itemss+median, nn-median);
        cilk_sync;
      }
      numItems = left->actualSize()+right->actualSize();
    }
  }

  void insert(objT** itemss, intT nn) {
    if(nn <= 0) return;
    counter += nn;
    if(nn < 2000) boundingBoxExpandSerial(itemss, nn);
    else boundingBoxExpandParallel(itemss, nn);

    if(isLeaf()) {//small
      if(numItems+nn >= capacity) resizeSerial((numItems+nn)*2);
      for(intT i=0; i<nn; ++i) {items[i+numItems]=itemss[i];}
      objT** newItems = &items[numItems];
      numItems += nn;
      if (actualSize() > leafSize) {
        floatT xM = -1;
        for (int kk=0; kk<dim; ++kk) {
          if (pMax[kk]-pMin[kk]>xM) {
            xM = pMax[kk]-pMin[kk];
            k = kk;}}
        xM = (pMax[k]+pMin[k])/2;
        split = splitterT(k, xM);
        compactSerial();
        intT median = splitItemSerial(items, numItems, split);
        left = new nodeT(items, median);
        right = new nodeT(items+median, numItems-median);
        makeInternal();
      }
    } else {//internal
      numItems += nn;
      intT median = splitItemSerial(itemss, nn, split);
      if(nn < 2000) {
        left->insert(itemss, median);
        right->insert(itemss+median, nn-median);
      } else {
        intT offset = left->size()+median;
        cilk_spawn left->insert(itemss, median);
        right->insert(itemss+median, nn-median);
        cilk_sync;
      }
    }
  }

  void printNodeRecurse(intT indent) {
    for(intT i=0; i<indent; ++i) cout << "- ";
    cout << "size " << actualSize() << "/" << size() << endl;
    if(!isLeaf()) {
      left->printNodeRecurse(indent+1);
      right->printNodeRecurse(indent+1);}
  }

  void rangeNeighbor(pointT pMin1, pointT pMax1, floatT r, vector<objT*>* accum) {
    int relation = boxCompare(pMin1, pMax1, pMin, pMax);
    if (relation == boxInclude) {
      if (isLeaf()) {
        for(intT i=0; i<numItems; ++i) if(items[i]) accum->push_back(items[i]);
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, accum);
        right->rangeNeighbor(pMin1, pMax1, r, accum);}
    } else if (relation == boxOverlap) {
      if (isLeaf()) {
        for(intT i=0; i<numItems; ++i) {
          if (items[i] && itemInBox(pMin1, pMax1, items[i])) accum->push_back(items[i]);}
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, accum);
        right->rangeNeighbor(pMin1, pMax1, r, accum);}
    }
  }

  template<class func, class func2>
  void rangeNeighbor(pointT pMin1, pointT pMax1, floatT r, func& f, func2& fStop) {
    if (fStop()) return;
    int relation = boxCompare(pMin1, pMax1, pMin, pMax);
    if (relation == boxInclude) {
      if (isLeaf()) {
        for(intT i=0; i<numItems; ++i) {
          if (items[i] && f(items[i])) break;
        }
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, f, fStop);
        right->rangeNeighbor(pMin1, pMax1, r, f, fStop);}
    } else if (relation == boxOverlap) {
      if (isLeaf()) {
        for(intT i=0; i<numItems; ++i) {
          if (items[i] && itemInBox(pMin1, pMax1, items[i])) {
            if (f(items[i])) break;}
        }
      } else {
        left->rangeNeighbor(pMin1, pMax1, r, f, fStop);
        right->rangeNeighbor(pMin1, pMax1, r, f, fStop);}
    }
  }
};

template<int dim, class objT>
class dynamicKdTree {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef dynamicKdNode<dim, objT> nodeT;
  nodeT *root;

  public:
  dynamicKdTree() {
    root=NULL;
  }
  dynamicKdTree(objT* P, intT n) {
    objT** items = newA(objT*, n);
    par_for(intT i=0; i<n; ++i) {items[i]=&P[i];}
    if (n>2000) {
      objT** scratch = newA(objT*, n);
      intT* flags = newA(intT, n);
      root = new nodeT(items, n, scratch, flags);
      free(scratch);
      free(flags);
    } else {
      root = new nodeT(items, n);}
    free(items);
  }
  void freeNodeRecurse(nodeT *node) {
    if(node->isLeaf()) {
      delete node;
    } else {
      freeNodeRecurse(node->leftChild());
      freeNodeRecurse(node->rightChild());
      delete node;}
  }
  ~dynamicKdTree() {
    freeNodeRecurse(root);
  }

  void erase(objT* PP, intT nn, bool rebuild=false) {
    if(!root) return;
    objT** PPP = newA(objT*, nn);
    par_for(intT i=0; i<nn; ++i) {PPP[i]=&PP[i];}
    root->erase(PPP, nn);
    free(PPP);
    if(rebuild) root->rebuildRecurse();
  }

  void insert(objT* PP, intT nn, bool rebuild=false) {
    if (!root) {
      objT** items = newA(objT*, nn);
      par_for(intT i=0; i<nn; ++i) {items[i]=&PP[i];}
      if (nn>2000) {
        objT** scratch = newA(objT*, nn);
        intT* flags = newA(intT, nn);
        root = new nodeT(items, nn, scratch, flags);
        free(scratch);
        free(flags);
      } else {
        root = new nodeT(items, nn);}
      free(items);
    } else {
      objT** PPP = newA(objT*, nn);
      par_for(intT i=0; i<nn; ++i) {PPP[i]=&PP[i];}
      root->insert(PPP, nn);
      free(PPP);
      if(rebuild) root->rebuildRecurse();
    }
  }

  vector<objT*>* rangeNeighbor(objT* query, floatT r) {
    vector<objT*>* accum = new vector<objT*>();
    if(!root) return accum;
    pointT pMin1 = pointT();
    pointT pMax1 = pointT();
    floatT* center = query->coordinate();
    for (int i=0; i<dim; ++i) {
      pMin1.updateX(i, center[i]-r);
      pMax1.updateX(i, center[i]+r);}
    root->rangeNeighbor(pMin1, pMax1, r, accum);
    return accum;
  }

  template<class func, class func2>
  void rangeNeighbor(objT* query, floatT r, func& f, func2& fStop) {
    if(!root) return;
    pointT pMin1 = pointT();
    pointT pMax1 = pointT();
    floatT* center = query->coordinate();
    for (int i=0; i<dim; ++i) {
      pMin1.updateX(i, center[i]-r);
      pMax1.updateX(i, center[i]+r);}
    root->rangeNeighbor(pMin1, pMax1, r, f, fStop);
  }

  intT size() {
    if(!root) return 0;
    return root->actualSize();
  }
  void printTree() {
    cout << "tree = " << endl;
    if(!root) cout << "empty" << endl;
    else root->printNodeRecurse(0);
  }
};

#endif
