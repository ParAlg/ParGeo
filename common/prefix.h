// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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

#ifndef PREFIX_H
#define PREFIX_H

#include "pbbs/utils.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"

template<class T, class F, class G>
void serial_prefix(T* A, intT n, F& process, G& cleanUp) {
  for(intT i=0; i<n; ++i) {
    if (process(A[i])) {
      cleanUp(A, i);}
  }
}

template<class T, class F, class G>
void parallel_prefix(T* A, intT n, F& process, G& cleanUp, bool verbose=false, intT* flag=NULL) {
  if (n < 2000) return serial_prefix(A, n, process, cleanUp);

  static const intT PREFIX = 1000;
  intT prefix;
  intT i = 0;
  intT conflict;

  bool freeFlag = false;
  if (!flag) {
    flag = newA(intT, n+1);
    freeFlag = true;}

  while (i < n) {
    if (i==0) prefix = PREFIX;//prefix < n is known
    else {
      prefix = min(i+i, n);//prefix to be taken as twice of i
      prefix = max(prefix, PREFIX);//prefix is too small, use larger, prefix < n is known
    }

    par_for (intT j=i; j<prefix; ++j) {
      if (process(A[j])) flag[j-i] = 1;//conflict
      else flag[j-i] = 0;
    }

    intT numBad = sequence::prefixSum(flag, 0, prefix-i);
    flag[prefix-i] = numBad;

    if (numBad > 0) {
      par_for(intT j=0; j<prefix-i; ++j) {
        if (flag[j]==0 && flag[j]!=flag[j+1]) conflict = j;}
      i += conflict;
      if (verbose) cout << "ci = " << i << endl;
      cleanUp(A, i);
    } else {
      i = prefix;
    }
  }

  if(freeFlag) free(flag);
}

#endif
