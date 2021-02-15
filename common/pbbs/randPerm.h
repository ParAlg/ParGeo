// This code is part of the paper "A Simple and Practical Linear-Work
// Parallel Algorithm for Connectivity" in Proceedings of the ACM
// Symposium on Parallelism in Algorithms and Architectures (SPAA),
// 2014.  Copyright (c) 2014 Julian Shun, Laxman Dhulipala and Guy
// Blelloch
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

#ifndef RAND_PERM_H
#define RAND_PERM_H

#include "parallel.h"
#include "utils.h"
#include "randPerm.h"
#include <iostream>
#include "sequence.h"
#include "speculative_for_2.h"
using namespace std;

inline unsigned int hashI(unsigned int a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

template <class E>
struct randPermStep {
  typedef pair<E,intT> pairInt;

  pairInt *dataCheck;
  intT* H;

  randPermStep(intT* _H, pairInt* _dataCheck) :
    H(_H), dataCheck(_dataCheck) {}

  inline bool reserve(intT i) {
    reserveLoc(dataCheck[i].second,i);
    reserveLoc(dataCheck[H[i]].second, i);
    return 1;
  }

  inline bool commit (intT i) {
    intT h = H[i];
    if(dataCheck[h].second == i) {
      if(dataCheck[i].second == i) {
	swap(dataCheck[i].first, dataCheck[h].first);
	dataCheck[h].second = intMax();
	return 1;
      }
      dataCheck[h].second = intMax();
    }
    return 0;
  }
};

template <class E>
void randPerm(E *A, intT n, intT r=50) {

  typedef pair<E,intT> pairInt;

  intT *H = newA(intT,n);
  //intT *check = newA(intT,n);
  pairInt* dataCheck = newA(pairInt,n);

  parallel_for (0, n, [&](intT i) {
			//H[i] = utils::hash(i)%(i+1);
			H[i] =i+hashI(i)%(n-i);
			dataCheck[i] = make_pair(A[i],i);
		      });

  randPermStep<E> rStep(H, dataCheck);
  speculative_for(rStep, 0, n, (r != -1) ? r : 50, 0);

  parallel_for (0, n, [&](intT i) {A[i] = dataCheck[i].first;});

  free(H); free(dataCheck);
}

#endif
