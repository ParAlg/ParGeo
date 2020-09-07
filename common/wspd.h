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

#ifndef WSPD
#define WSPD

#include "kdTree.h"
#include "geometry.h"

template<class nodeT, class opT>
inline void findPairSerial(nodeT *u, nodeT *v, opT* f) {
  if (!f->moveon(u, v)) return;
  if (u->wellSeparated(v)) {
    f->run(u, v);
  } else {
    if (u->isLeaf() && v->isLeaf()) {
      cout << "error, leaves not well separated, abort" << endl;
      abort();
    } else if (u->isLeaf()) {
      findPairSerial(v->L(), u, f);
      findPairSerial(v->R(), u, f);
    } else if (v->isLeaf()) {
      findPairSerial(u->L(), v, f);
      findPairSerial(u->R(), v, f);
    } else {
      if (u->lMax() > v->lMax()) {
        findPairSerial(u->L(), v, f);
        findPairSerial(u->R(), v, f);
      } else {
        findPairSerial(v->L(), u, f);
        findPairSerial(v->R(), u, f);}
    }}
}

template<class nodeT, class opT>
inline void findPairParallel(nodeT *u, nodeT *v, opT* f) {
  if (!f->moveon(u, v)) return;
  if (u->size() + v->size() < 2000) return findPairSerial(u,v,f);

  if (u->wellSeparated(v)) {
    f->run(u, v);//need to be thread safe
  } else {
    if (u->isLeaf() && v->isLeaf()) {
      cout << "error, leaves not well separated, abort" << endl;
      abort();
    } else if (u->isLeaf()) {
      cilk_spawn findPairSerial(v->L(), u, f);
      findPairSerial(v->R(), u, f);
      cilk_sync;
    } else if (v->isLeaf()) {
      cilk_spawn findPairSerial(u->L(), v, f);
      findPairSerial(u->R(), v, f);
      cilk_sync;
    } else {
      if (u->lMax() > v->lMax()) {
        cilk_spawn findPairSerial(u->L(), v, f);
        findPairSerial(u->R(), v, f);
        cilk_sync;
      } else {
        cilk_spawn findPairSerial(v->L(), u, f);
        findPairSerial(v->R(), u, f);
        cilk_sync;
      }
    }}
}

template<class nodeT, class opT>
inline void wspdSerial(nodeT *nd, opT *f) {
  if (!nd->isLeaf() && f->start(nd)) {
    wspdSerial(nd->L(), f);
    wspdSerial(nd->R(), f);
    findPairSerial<nodeT, opT>(nd->L(), nd->R(), f);}
}

template<class nodeT, class opT>
inline void wspdParallel(nodeT *nd, opT *f) {
  if (nd->size() < 2000) {
    wspdSerial(nd, f);
  } else if (!nd->isLeaf() && f->start(nd)) {
    cilk_spawn wspdParallel(nd->L(), f);
    wspdParallel(nd->R(), f);
    cilk_sync;
    findPairParallel<nodeT, opT>(nd->L(), nd->R(), f);}
}

#endif
