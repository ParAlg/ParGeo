// This code is part of the project "A Parallel Batch-Dynamic Data Structure
// for the Closest Pair Problem"
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

#ifndef SERIAL_HEAP_H
#define SERIAL_HEAP_H

/* Heap specifically implemented for closest pair <dist,u,v> maintenance.
   It keeps both a standard pq keyed on dist
   and an additional hashtable to retrieve dist just based on u.
   So it supports both get-min based on dist and erase/find based on u.
 */

#include "../common/pbbs/ndHash.h"
#include "../common/shared.h"

template<int dim> struct pointer;
template<int dim> struct pointerHash;
template<int dim> struct pointPairHeapCmp;

template<int dim>
struct serialHeap {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointer<dim> pointerT;
  typedef pointerHash<dim> hashT;
  typedef Table<hashT,intT> tableT;
  typedef pointPair<dim> pointPairT;
  typedef hashFloatToCell<dim> cellHashT;
  typedef set<pointPairT, pointPairHeapCmp<dim>> heapT;
  cellHashT* cellHash=NULL;
  tableT* table=NULL;
  heapT* heap=NULL;
  pointerT* buffer=NULL;
  intT bc=0;
  intT heapMax;

  serialHeap(intT heapMaxx, point<dim> pMin, floatT r): heapMax(heapMaxx) {
    cellHash = new cellHashT(pMin, r);
    table = new tableT(heapMax*2, hashT(cellHash));
    buffer = newA(pointerT, heapMax);//tune heap max
    for(intT i=0; i<heapMax; ++i) buffer[i]=pointerT();
    heap = new heapT();
  }

  ~serialHeap() {
    if(cellHash) delete cellHash;
    if(table) table->del();
    if(heap) delete heap;
    if(buffer) free(buffer);}

  inline void doubleBuffer() {
    auto newBuffer = newA(pointerT, heapMax*2);
    bc=0;
    table->clear();
    for(intT i=0; i<heapMax; ++i) {
      if(!buffer[i].isEmpty()) {
        newBuffer[bc]=buffer[i];
        table->insert(&newBuffer[bc]);
        bc++;}
    }
    for(intT i=bc; i<heapMax*2; ++i) {newBuffer[i]=pointerT();}
    heapMax*=2;
    swap(buffer, newBuffer);
    free(newBuffer);
  }

  inline pointerT* writeBuffer(pointerT pp) {
    if (bc>heapMax-1) doubleBuffer();
    buffer[bc] = pp;
    return &buffer[bc++];}

  //user's responsibility to ensure no repeats, todo, handle duplication
  void insert(pointPairT pp) {
    auto tmp = writeBuffer(pointerT(pp.u, pp.dist));
    table->insert(tmp);
    heap->insert(pp);
  }

  bool erase(pointT u) {
    auto u0 = table->find(&u);
    if (u0->isEmpty()) {
      return false;
    } else {
      table->deleteVal(&u);
      heap->erase(pointPairT(u,u,u0->dist));
      u0->setEmpty();
      return true;}
  }

  //if pp.dist < heap(pp.u).dist, decrease key
  bool decreaseKey(pointPairT pp) {
    auto u0 = table->find(&pp.u);
    if (u0->isEmpty()) {
      return false;
    } else if (pp.dist<u0->dist) {
      heap->erase(pointPairT(pp.u,pp.u,u0->dist));
      heap->insert(pp);
      u0->dist = pp.dist;
      return true;}
    return false;
  }

  //replace heap(pp.u) with pp
  bool replace(pointPairT pp) {
    auto u0 = table->find(&pp.u);
    if (u0->isEmpty()) {
      return false;
    } else {
      heap->erase(pointPairT(pp.u,pp.u,u0->dist));
      heap->insert(pp);
      return true;}
  }

  floatT dist(pointT u) {
    auto u0 = table->find(&u);
    if (u0->isEmpty()) return pointerT::floatMax;
    else return u0->dist;}

  intT size() {return heap->size();}
  bool consistent() {return heap->size()==table->count();}
  pointPairT getMin() {return *(heap->begin());}
};

// *************************************************************
//    Pointer <u,dist>
// *************************************************************

//points from u to distance
template<int dim>
struct pointer {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointer<dim> pointerT;
  static constexpr double floatMax = numeric_limits<floatT>::max();
  pointT u;
  floatT dist;
  pointer(pointT uu, floatT distt):u(uu),dist(distt) {}
  pointer():dist(floatMax) {u=pointT();}
  bool isEmpty() {return dist>=floatMax-1 && u.isEmpty();}
  bool setEmpty() {return dist=floatMax; u=pointT();}
  floatT* coordinate() {return u.x;}
  bool operator==(pointerT p2) {
    return samePoint<dim>(u.x, p2.u.x);
  }
};

template<int dim>
struct pointerHash {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef hashFloatToCell<dim> hashFunc;
  typedef pointer<dim> pointerT;
  typedef pointerT* eType;
  typedef pointT* kType;
  hashFunc* hashF;
  eType e;
  pointerHash(hashFunc* hashFF):hashF(hashFF) {
    e = new pointerT();
  }
  eType empty() {return e;}
  kType getKey(eType v) {return &v->u;}
  uintT hash(kType c) {return hashF->hash(c->coordinate());}
  int cmp(kType c1, kType c2) {
    if (c1->isEmpty() || c2->isEmpty()) return 1;
    return hashF->comparePoint(c1->coordinate(), c2->coordinate());}
  bool replaceQ(eType c1, eType c2) {return 1;}
};

template<int dim>
struct pointPairHeapCmp {
  typedef pointPair<dim> pointPairT;
  bool operator() (pointPairT i, pointPairT j){
    return i<j;
  }
};

#endif
