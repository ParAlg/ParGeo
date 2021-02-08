// This code is part of the Pargeo Library
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

#include <algorithm>
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "quick.h"
#include "graham.h"
#include "line.h"
#include "randHull.h"
using namespace std;
using namespace sequence;

_seq<intT> hull(point2d* P, intT n) {
  //0: quickHull, 1: NA, 2: graham, 3: lineSweep, 4: randHull
  static const intT baseMethod = 4;
  auto procs = getWorkers();
  timing t; t.start();
  timing t1; t1.start();

  if (n < 2000) {
    _seq<intT> CH;
    CH.A = newA(intT, n);
    if (baseMethod == 0) { CH.n = quickHullSerial(P, n, CH.A);}
    else if (baseMethod == 1) { abort(); }
    else if (baseMethod == 2) { CH.n = grahamScanSerial(P, n, CH.A);}
    else if (baseMethod == 3) { CH.n = lineSweepSerial(P, n, CH.A);}
    else if (baseMethod == 4) { CH.n = randHullSerial(P, n, CH.A);}
    else {
      cout << "Error, wrong method number, abort." << endl; abort();}
#ifndef SILENT
    cout << "serial-hull-time = " << t.stop() << endl;
#endif
    check(P, n, CH.A, CH.n);
    return CH;
  }

  intT blk = n/procs;
  intT M[procs+1];

  intT* I = newA(intT, n);
  parallel_for(0, n, [&](intT i) {I[i] = i;});
  intT* IOut = newA(intT, n);

  parallel_for(0, procs,
	       [&](intT i) {
		 intT o = blk*i;
		 if (i != procs-1) {
		   if (baseMethod == 0) M[i] = quickHullSerial(P, blk, I+o, true);
		   else if (baseMethod == 1) abort();
		   else if (baseMethod == 2) M[i] = grahamScanExternalSerial(P, I+o, blk, IOut+o);
		   else if (baseMethod == 3) M[i] = lineSweepExternalSerial(P, I+o, blk, IOut+o);
		   else if (baseMethod == 4) {
		     IOut[o+0] = I[o+0];
		     IOut[o+1] = I[o+1];
		     M[i] = randHullExternalSerial(P, I+o, blk, IOut+o, 2);
		   } else {
		     cout << "Error, wrong method number, abort." << endl; abort();}
		 } else {
		   intT blkLast = max(blk, n-blk*(procs-1));
		   if (baseMethod == 0) M[i] = quickHullSerial(P, blkLast, I+o, true);
		   else if (baseMethod == 1) abort();
		   else if (baseMethod == 2) M[i] = grahamScanExternalSerial(P, I+o, blkLast, IOut+o);
		   else if (baseMethod == 3) M[i] = lineSweepExternalSerial(P, I+o, blkLast, IOut+o);
		   else if (baseMethod == 4) {
		     IOut[o+0] = I[o+0];
		     IOut[o+1] = I[o+1];
		     M[i] = randHullExternalSerial(P, I+o, blkLast, IOut+o, 2);
		   } else {
		     cout << "Error, wrong method number, abort." << endl; abort();}
		 }
	       }, 1);
  if (baseMethod == 0 || baseMethod == 4) swap(I, IOut);
#ifndef SILENT
  cout << "subhull-time = " << t1.next() << endl;
#endif

  //Prefix-sum M
  for(intT i=1; i<procs+1; ++i) {
    M[i] += M[i-1];}
  for(intT i=procs; i>0; --i) {
    M[i] = M[i-1];}
  M[0] = 0;
#ifndef SILENT
  cout << "subhull-size = " << M[procs] << endl;
#endif

  parallel_for(0, procs, [&](intT i) {
			   intT o = blk*i;
			   for(intT j=M[i]; j<M[i+1]; ++j) {
			     I[j] = IOut[o+j-M[i]];
			   }
			 }, 1);

  intT m2;
  if (M[procs] < 2000) {
    if (baseMethod == 0 || baseMethod == 4) {//todo
      m2 = quickHullSerial(P, M[procs], I, true);
      swap(I, IOut);
    } else if (baseMethod == 1) { abort(); }
    else if (baseMethod == 2) {
      m2 = grahamScanExternalSerial(P, I, M[procs], IOut);
    } else if (baseMethod == 3) {
      m2 = lineSweepExternalSerial(P, I, M[procs], IOut);
    } else if (baseMethod == 4) {
      //todo bug here, fine with other method
      IOut[0] = I[0];
      IOut[1] = I[1];
      m2 = randHullExternalSerial(P, I, M[procs], IOut, 2);
      swap(I, IOut);
    } else {
      cout << "Error, wrong method number, abort." << endl; abort();}
  } else {
    m2 = quickHullParallel(P, M[procs], I, IOut, true);
    swap(I, IOut);
  }

#ifndef SILENT
  cout << "merge-time = " << t1.next() << endl;
  cout << "total-hull-time = " << t.stop() << endl;
  check(P, n, IOut, m2);
#else
  cout << t.stop() << endl;
#endif

  free(I);
  return _seq<intT>(IOut, m2);
}
