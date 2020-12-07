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

template <class nodeT>
struct wsp {
  nodeT* u;
  nodeT* v;
  wsp(nodeT* uu, nodeT* vv): u(uu), v(vv) {}
};

template<class nodeT, class opT>
inline void findPairSerial(nodeT *u, nodeT *v, opT* f, floatT s) {
  if (!f->moveon(u, v)) return;

  if (f->wellSeparated(u, v, s)) {
    f->run(u, v);
  } else {
    if (u->isLeaf() && v->isLeaf()) {
      cout << "error, leaves not well separated, abort" << endl;
      abort();
    } else if (u->isLeaf()) {
      findPairSerial(v->L(), u, f, s);
      findPairSerial(v->R(), u, f, s);
    } else if (v->isLeaf()) {
      findPairSerial(u->L(), v, f, s);
      findPairSerial(u->R(), v, f, s);
    } else {
      if (u->lMax() > v->lMax()) {
        findPairSerial(u->L(), v, f, s);
        findPairSerial(u->R(), v, f, s);
      } else {
        findPairSerial(v->L(), u, f, s);
        findPairSerial(v->R(), u, f, s);}
    }}
}

template<class nodeT, class opT>
inline void findPairParallel(nodeT *u, nodeT *v, opT* f, floatT s) {
  if (!f->moveon(u, v)) return;
  if (u->size() + v->size() < 2000) return findPairSerial(u,v,f,s);

  if (f->wellSeparated(u, v, s)) {
    f->run(u, v);//need to be thread safe
  } else {
    if (u->isLeaf() && v->isLeaf()) {
      cout << "error, leaves not well separated, abort" << endl;
      abort();
    } else if (u->isLeaf()) {
      par_do([&](){findPairParallel(v->L(), u, f, s);},
	     [&](){findPairParallel(v->R(), u, f, s);});
    } else if (v->isLeaf()) {
      par_do([&](){findPairParallel(u->L(), v, f, s);},
	     [&](){findPairParallel(u->R(), v, f, s);});
    } else {
      if (u->lMax() > v->lMax()) {
	par_do([&](){findPairParallel(u->L(), v, f, s);},
	       [&](){findPairParallel(u->R(), v, f, s);});
      } else {
	par_do([&](){findPairParallel(v->L(), u, f, s);},
	       [&](){findPairParallel(v->R(), u, f, s);});
      }
    }}
}

template<class nodeT, class opT>
inline void wspdSerial(nodeT *nd, opT *f, floatT s=2) {
  if (!nd->isLeaf() && f->start(nd)) {
    wspdSerial(nd->L(), f, s);
    wspdSerial(nd->R(), f, s);
    findPairSerial<nodeT, opT>(nd->L(), nd->R(), f, s);}
}

template<class nodeT, class opT>
inline void wspdParallel(nodeT *nd, opT *f, floatT s=2) {
  if (nd->size() < 2000) {
    wspdSerial(nd, f, s);
  } else if (!nd->isLeaf() && f->start(nd)) {
    par_do([&](){wspdParallel(nd->L(), f, s);},
	   [&](){wspdParallel(nd->R(), f, s);});
    findPairParallel<nodeT, opT>(nd->L(), nd->R(), f, s);}
}

#endif
