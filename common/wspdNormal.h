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

#ifndef WSPD_NORMAL
#define WSPD_NORMAL

#include <vector>
#include "parBuf.h"
#include "wspd.h"

template <class nodeT>
struct wspdNormalSerial {
  typedef wsp<nodeT> pType;
  vector<pType> *out;

  wspdNormalSerial(vector<pType> *outt) : out(outt) {}

  inline void run(nodeT *u, nodeT *v) {out->emplace_back(u, v);}
  inline bool moveon(nodeT *u, nodeT *v) {return true;}
  inline bool start(nodeT *u) { return true; }
  inline bool wellSeparated(nodeT *u, nodeT *v, int s) {return u->wellSeparated(v, s);}
};

template <class nodeT>
struct wspdNormalParallel {
  typedef wsp<nodeT> pType;
  typedef parBuf<pType> bufT;
  bufT **out;

  wspdNormalParallel(intT n) {
    int P = getWorkers();
    out = newA(bufT*, P);
    par_for(int p=0; p<P; ++p) {
      out[p] = new bufT(n/P);
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
  inline bool start(nodeT *t_u) { return true; }
  inline bool wellSeparated(nodeT *u, nodeT *v, int s) {return u->wellSeparated(v, s);}
};

#endif
