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
#include "pbbs/gettime.h"

template<class T, class F, class G>
void serial_prefix(T* A, intT n, F& process, G& cleanUp) {
  for(intT i=0; i<n; ++i) {
    if (process(A[i])) {
      cleanUp(A, i);}
  }
}

template<class T, class F, class G>
void parallel_prefix(T* A, intT n, F& process, G& cleanUp, bool verbose=false, intT* flag=NULL) {
  static const intT PREFIX = 200000;
  static const intT THRESH = 500000;
  if (n < THRESH) return serial_prefix(A, n, process, cleanUp);

  intT prefix;
  intT i = 0;
  intT conflict;

  bool freeFlag = false;
  if (!flag) {
    flag = newA(intT, n+1);
    freeFlag = true;}

  timing t0;
  intT serCount = 0;
  intT serSize = 0;
  floatT serTime = 0;
  intT parCount = 0;
  intT parSize = 0;
  floatT parTime = 0;
  floatT parTime1 = 0;
  floatT parTime2 = 0;
  floatT parTime3 = 0;

  while (i < n) {
    if (i==0) prefix = PREFIX;//prefix < n is known
    else prefix = min(i*3, n);//prefix to be taken as 3*i

    //if (verbose) cout << "prefix = " << prefix << endl;

    if (prefix-i < THRESH) {
      if(verbose) t0.start();
      serCount += 1;
      serSize += prefix - i;
      while (i < prefix) {
	if (process(A[i])) cleanUp(A, i);
	i ++;
      }
      if(verbose) serTime += t0.stop();
    } else {
      if(verbose) t0.start();
      parCount += 1;
      parSize += prefix - i;

      par_for (intT j=i; j<prefix; ++j) {
	if (process(A[j])) flag[j-i] = 1;//conflict
	else flag[j-i] = 0;
      }
      if(verbose) parTime += t0.next();

      intT numBad = sequence::prefixSum(flag, 0, prefix-i);
      flag[prefix-i] = numBad;
      if(verbose) parTime1 += t0.stop();

      if (numBad > 0) {
	if(verbose) t0.start();
	par_for(intT j=0; j<prefix-i; ++j) {
	  if (flag[j]==0 && flag[j]!=flag[j+1]) conflict = j;}
	i += conflict;
	if(verbose) parTime2 += t0.next();
	//if (verbose) cout << "ci = " << i << endl;
	cleanUp(A, i);
	if(verbose) parTime3 += t0.stop();
      } else {
	i = prefix;
      }

    }
  }

  if(verbose) {
    cout << "serial prefix stats" << endl;
    cout << " count = " << serCount << endl;
    cout << " avg size = " << serSize/(floatT)serCount << endl;
    cout << " total time = " << serTime << endl;
    cout << "parallel prefix stats" << endl;
    cout << " count = " << parCount << endl;
    cout << " avg size = " << parSize/(floatT)parCount << endl;
    cout << " total time 1 = " << parTime << endl;
    cout << " total time 2 = " << parTime1 << endl;
    cout << " total time 3 = " << parTime2 << endl;
    cout << " total time 4 = " << parTime3 << endl;
  }

  if(freeFlag) free(flag);
}

#endif
