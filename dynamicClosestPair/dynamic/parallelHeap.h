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

#ifndef PARALLEL_HEAP_H
#define PARALLEL_HEAP_H

#include "serialHeap.h"
#include "PAM/c++/pam.h"
#include "pbbslib/seq.h"
#include "pbbslib/parallel.h"
#include "pbbslib/get_time.h"

// *************************************************************
//    PAM connector, rb-tree + map
// *************************************************************

//pam tree node for type pp (e.g. pointPair<dim>)
template<class pp>
struct nodePP {
  using key_t = pp;
  using aug_t = double;
  using entry_t = pp;
  using val_t = entry_t;
  static inline key_t get_key(const entry_t& e) {return e;}
  static inline val_t get_val(const entry_t& e) {return e;}
  static inline void set_val(entry_t& e, const val_t& v) {
    e.u = v.u;
    e.v = v.v;
    e.dist = v.dist;}
  static inline aug_t from_entry(const entry_t& e) {//unused
    return e.dist;}
  static inline bool comp(key_t a, key_t b) { return a<b;}
  static aug_t get_empty() { return entry_t::floatMax;}//unused
  static aug_t combine(aug_t a, aug_t b) { return b;}//unused
};

template <class _Entry, class Balance=weight_balanced_tree>
using aug_map_pp =
  aug_map_<nodePP<_Entry>,
  typename Balance::template
  balance<aug_node<typename Balance::data,
		   nodePP<_Entry>>>>;

// *************************************************************
//    parallel heap wrapper
// *************************************************************

template<int dim> struct pointer;
template<int dim> struct pointerHash;

//todo serial basecases
template<int dim>
struct parallelHeap {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointer<dim> pointerT;
  typedef pointerHash<dim> hashT;
  typedef Table<hashT,intT> tableT;
  typedef pointPair<dim> pointPairT;
  typedef hashFloatToCell<dim> cellHashT;
  typedef aug_map_pp<pointPair<dim>> heapT;
  typedef set<pointPairT, pointPairHeapCmp<dim>> stlSet;
  static const bool useSerial = false;
  cellHashT* cellHash=NULL;
  tableT* table=NULL;
  heapT heap;
  stlSet heap2;
  pointerT* buffer=NULL;
  intT bc=0;
  intT heapMax;

  struct statistics {
    double insertBufferTime = 0;
    double insertTableTime = 0;
    double insertHeapTime = 0;
    double insertTotal = 0;
    double insertSerialTotal = 0;
    intT insertCount = 0;
    intT insertSerialCount = 0;
    double eraseBufferTableTime = 0;
    double eraseHeapTime = 0;
    double eraseTotal = 0;
    double eraseSerialTotal = 0;
    intT eraseCount = 0;
    intT eraseSerialCount = 0;
  };
  struct statistics stats;
  void printStats() {
    std::cout << std::setprecision(3);
    cout << ">>> heap stats " << endl << endl;
    cout << "insert parallel count = " << stats.insertCount << endl;
    cout << "insert parallel time total = " << stats.insertTotal << endl;
    cout << "  buffer time = " << stats.insertBufferTime;
    cout << " (" << intT(100*stats.insertBufferTime/stats.insertTotal) << "%)" << endl;
    cout << "  table time = " << stats.insertTableTime;
    cout << " (" << intT(100*stats.insertTableTime/stats.insertTotal) << "%)" << endl;
    cout << "  heap time = " << stats.insertHeapTime;
    cout << " (" << intT(100*stats.insertHeapTime/stats.insertTotal) << "%)" << endl;
    cout << "insert serial count = " << stats.insertSerialCount << endl;
    cout << "insert serial time total = " << stats.insertSerialTotal << endl;
    cout << endl;
    cout << "erase parallel count = " << stats.eraseCount << endl;
    cout << "erase parallel time total = " << stats.eraseTotal << endl;
    cout << "  buffer+table time = " << stats.eraseBufferTableTime;
    cout << " (" << intT(100*stats.eraseBufferTableTime/stats.eraseTotal) << "%)" << endl;
    cout << "  heap time = " << stats.eraseHeapTime;
    cout << " (" << intT(100*stats.eraseHeapTime/stats.eraseTotal) << "%)" << endl;
    cout << "erase serial count = " << stats.eraseSerialCount << endl;
    cout << "erase serial time total = " << stats.eraseSerialTotal << endl;
    cout << endl;
    std::cout << std::setprecision(6);
  }

  parallelHeap(intT heapMaxx, point<dim> pMin, floatT r): heapMax(heapMaxx) {
    cellHash = new cellHashT(pMin, r);
    table = new tableT(heapMax*2, hashT(cellHash));
    buffer = newA(pointerT, heapMax);//tune heap max
    for(intT i=0; i<heapMax; ++i) buffer[i]=pointerT();
    if(useSerial) heap2 = stlSet();
    else heap = heapT();
  }

  ~parallelHeap() {
    if(cellHash) delete cellHash;
    if(table) table->del();
    if(buffer) free(buffer);}

  //accommodates at least nn additional items
  inline void doubleBufferSerial(intT nn) {
    intT newSize = (nn+bc)*2;
    auto newBuffer = newA(pointerT, newSize);
    intT ii=0;
    table->clear();
    for(intT i=0; i<bc; ++i) {
      if(!buffer[i].isEmpty()) {
        newBuffer[ii] = buffer[i];
        table->insert(&newBuffer[ii]);
        ii++;}
    }
    for(intT i=ii; i<newSize; ++i) newBuffer[i]=pointerT();
    heapMax = newSize;
    swap(buffer, newBuffer);
    free(newBuffer);
  }

  inline void doubleBuffer(intT nn) {
    intT* flag = newA(intT, bc+1);
    auto fMark = [&](long i){
                   if(!buffer[i].isEmpty()) flag[i] = 1;
                   else flag[i] = 0;};
    parallel_for(0, bc, fMark, 0);

    intT nElem = sequence::prefixSum(flag, 0, bc);
    flag[bc] = nElem;
    intT newSize = (nn+nElem)*2;
    auto newBuffer = newA(pointerT, newSize);
    table->clear();
    auto fInsert = [&](long i){
                     if(flag[i] != flag[i+1]) {
                       newBuffer[flag[i]] = buffer[i];
                       table->insert(&newBuffer[flag[i]]);}
                   };
    parallel_for(0, bc, fInsert, 0);
    bc = nElem;
    parallel_for(bc, newSize, [&](long i){newBuffer[i]=pointerT();}, 0);
    heapMax = newSize;
    swap(buffer, newBuffer);
    free(newBuffer);
    free(flag);
  }

  inline pointerT* writeBufferSerial(pointPairT* PP, intT nn) {
    if (bc+nn>heapMax-1) doubleBufferSerial(nn);
    for(intT i=0; i<nn; ++i){
      buffer[bc+i] = pointerT(PP[i].u, PP[i].dist);}
    auto writeStart = &buffer[bc];
    bc += nn;
    return writeStart;
  }

  inline pointerT* writeBuffer(pointPairT* PP, intT nn) {
    if (bc+nn>heapMax-1) doubleBuffer(nn);
    parallel_for(0, nn, [&](long i){
                          buffer[bc+i] = pointerT(PP[i].u, PP[i].dist);}, 0);
    auto writeStart = &buffer[bc];
    bc += nn;
    return writeStart;
  }

  //user's responsibility to handle duplicates, todo, handle duplicates
  void insertSerial(pointPairT* PP, intT nn) {
    stats.insertSerialCount ++;
    timer tt;
    tt.start();
    auto tmp = writeBufferSerial(PP, nn);
    for(intT i=0; i<nn; ++i) {
      table->insert(&tmp[i]);
      if(useSerial) heap2.insert(PP[i]);
      else heap = heapT::insert(move(heap), PP[i]);
    }
    stats.insertSerialTotal += tt.stop();
  }

  //mutates PP
  void insert(pointPairT* PP, intT nn) {
    if(nn<400) return insertSerial(PP, nn);

    stats.insertCount ++;
    timer t0, tt;
    tt.start();

    t0.start();
    auto tmp = writeBuffer(PP, nn);
    stats.insertBufferTime += t0.stop();

    t0.start();
    auto pInsert = [&](long i) {
                     table->insert(&tmp[i]);};
    parallel_for(0, nn, pInsert, 0);
    stats.insertTableTime += t0.stop();

    t0.start();
    if(useSerial) {
      for(intT i=0; i<nn; ++i) heap2.insert(PP[i]);
    } else {
      auto PP2 = pbbs::sequence<pointPairT>(nn);
      parallel_for(0, nn, [&](long i){PP2[i] = PP[i];}, 0);
      heap = heapT::multi_insert(move(heap), PP2);
    }
    stats.insertHeapTime += t0.stop();
    stats.insertTotal += tt.stop();
  }

  void eraseSerial(pointT* P, intT nn) {
    stats.eraseSerialCount ++;
    timer tt;
    tt.start();
    for (intT i=0; i<nn; ++i) {
      auto u0 = table->find(&P[i]);
      if (!u0->isEmpty()) {
        table->deleteVal(&P[i]);
        auto tmp = pointPairT(P[i],P[i],u0->dist);
        if(useSerial) heap2.erase(tmp);
        else heap = heapT::remove(move(heap), tmp);
        u0->setEmpty();}
    }
    stats.eraseSerialTotal += tt.stop();
  }

  void erase(pointT* P, intT nn) {
    if (nn < 400) return eraseSerial(P, nn);

    stats.eraseCount ++;
    timer t0, tt;
    tt.start();

    t0.start();
    auto flag = newA(intT, nn+1);
    auto fMark = [&](long i) {
                   auto u0 = table->find(&P[i]);
                   if (u0->isEmpty()) flag[i] = 0;
                   else {
                     flag[i] = 1;}
                 };
    parallel_for(0, nn, fMark, 0);
    intT nElem = sequence::prefixSum(flag, 0, nn);
    flag[nn] = nElem;
    auto tmp = pbbs::sequence<pointPairT>(nElem);
    auto fTableDel = [&](long i){
                       if (flag[i]!=flag[i+1]) {
                         auto u = P[i];
                         auto u0 = table->find(&u);
                         table->deleteVal(&u);
                         tmp[flag[i]] = pointPairT(u,u,u0->dist);
                         u0->setEmpty();}
                     };
    parallel_for(0, nn, fTableDel, 0);
    stats.eraseBufferTableTime += t0.stop();

    t0.start();
    if(useSerial) {
      for(intT i=0; i<nElem; ++i) heap2.erase(tmp[i]);
    } else {
      heap = heapT::multi_delete(move(heap), tmp);
    }
    stats.eraseHeapTime += t0.stop();
    free(flag);

    stats.eraseTotal += tt.stop();
  }

  floatT dist(pointT u) {
    auto u0 = table->find(&u);
    if (u0->isEmpty()) return pointerT::floatMax;
    else return u0->dist;}

  intT size() {
    if(useSerial) return heap2.size();
    else return heap.size();
  }
  bool consistent() {return heap.size()==table->count();}
  pointPairT getMin() {
    if(useSerial) return *(heap2.begin());
    else return *heap.select(0);
  }
};

#endif
