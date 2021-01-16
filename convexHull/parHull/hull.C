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
#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "quick.h"
#include "gift.h"
#include "graham.h"
#include "line.h"
#include "randHull.h"
using namespace std;
using namespace sequence;

_seq<intT> hull(point2d* P, intT n) {
  static const bool sortPoint=false;
  //0: quickHull, 1: giftWrap, 2: graham, 3: lineSweep, 4: randHull
  static const intT baseMethod = 0;
  auto procs = getWorkers();
  timing t; t.start();
  timing t1; t1.start();

  if (n < 2000) {
    _seq<intT> CH;
    if (baseMethod == 0) CH = quickHullSerial(P, n);
    else if (baseMethod == 1) CH = giftWrapSerial(P, n);
    else if (baseMethod == 2) CH = grahamScanSerial(P, n);
    else if (baseMethod == 3) CH = lineSweepSerial(P, n);
    else if (baseMethod == 4) CH = randHullSerial(P, n);
    else {
      cout << "Error, wrong method number, abort." << endl; abort();}
#ifndef SILENT
    cout << "serial-hull-time = " << t.stop() << endl;
#endif
    check(P, n, CH.A, CH.n);
    return CH;
  }

  if (sortPoint) {
    auto xLess = [&](point2d p1, point2d p2) {
		   return p1.x()<p2.x();};
    sampleSort(P, n, xLess);
#ifndef SILENT
    cout << "sort-time = " << t1.next() << endl;
#endif
  }

  intT blk = n/procs;
  intT M[procs+1];

  intT* I = newA(intT, n);

  parallel_for(0, procs,
	       [&](intT i) {
		 intT o = blk*i;
		 if (i != procs-1) {
		   if (baseMethod == 0) M[i] = quickHullSerial(P+o, blk, I+o).n;
		   else if (baseMethod == 1) M[i] = giftWrapSerial(P+o, blk, I+o).n;
		   else if (baseMethod == 2) M[i] = grahamScanSerial(P+o, blk, I+o).n;
		   else if (baseMethod == 3) M[i] = lineSweepSerial(P+o, blk, I+o).n;
		   else if (baseMethod == 4) M[i] = randHullSerial(P+o, blk, I+o).n;
		   else {
		     cout << "Error, wrong method number, abort." << endl; abort();}
		   //check(P+o, blk, I+o, M[i]);
		 } else {
		   intT blkLast = max(blk, n-blk*(procs-1));
		   if (baseMethod == 0) M[i] = quickHullSerial(P+o, blkLast, I+o).n;
		   else if (baseMethod == 1) M[i] = giftWrapSerial(P+o, blkLast, I+o).n;
		   else if (baseMethod == 2) M[i] = grahamScanSerial(P+o, blkLast, I+o).n;
		   else if (baseMethod == 3) M[i] = lineSweepSerial(P+o, blkLast, I+o).n;
		   else if (baseMethod == 4) M[i] = randHullSerial(P+o, blkLast, I+o).n;
		   else {
		     cout << "Error, wrong method number, abort." << endl; abort();}
		   //check(P+o, blkLast, I+o, M[i]);
		 }
	       }, 1);
#ifndef SILENT
  cout << "subhull-time = " << t1.next() << endl;
#endif
  for(intT i=1; i<procs; ++i) {
    M[i] += M[i-1];}
  for(intT i=procs; i>0; --i) {
    M[i] = M[i-1];}
  M[0] = 0;
#ifndef SILENT
  cout << "subhull-size = " << M[procs] << endl;
#endif
  point2d* PP = newA(point2d, M[procs]);

  parallel_for(0, procs, [&](intT i) {
			   intT o = blk*i;
			   for(intT j=M[i]; j<M[i+1]; ++j) {
			     PP[j] = P[o+I[o+(j-M[i])]];}
			 }, 1);

  intT m2;
  if (M[procs] < 2000) {
    if (baseMethod == 0) {
      m2 = quickHullSerial(PP, M[procs], I).n;
    } else if (baseMethod == 1) {
      m2 = giftWrapSerial(PP, M[procs], I).n;
    } else if (baseMethod == 2) {
      m2 = grahamScanSerial(PP, M[procs], I).n;
    } else if (baseMethod == 3) {
      m2 = lineSweepSerial(PP, M[procs], I).n;
    } else if (baseMethod == 4) {
      m2 = randHullSerial(PP, M[procs], I).n;
    } else {
      cout << "Error, wrong method number, abort." << endl; abort();}
  } else {
    _seq<intT> CH;
    if (baseMethod == 0 || baseMethod == 2 || baseMethod == 3 || baseMethod == 4) {
      CH = quickHullParallel(PP, M[procs]);
    } else if (baseMethod == 1) {
      CH = giftWrapParallel(PP, M[procs]);
    } else {
      cout << "Error, wrong method number, abort." << endl; abort();}
    m2 = CH.n;
    free(I);
    I = CH.A;
  }
#ifndef SILENT
  cout << "merge-time = " << t1.next() << endl;
  cout << "total-hull-time = " << t.stop() << endl;
  //hull is in (PP, I), size is m2
  check(P, n, I, m2, PP);
#else
  cout << t.stop() << endl;
#endif

  free(I);
  free(PP);
  return _seq<intT>();
}
