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
#include "serialQuick.h"
#include "parQuick.h"
using namespace std;
using namespace sequence;

_seq<intT> hull(point2d* P, intT n) {
  static const bool sortPoint=false;
  auto procs = getWorkers();
  timing t; t.start();
  timing t1; t1.start();

  if (n < 2000) {
    intT* I = newA(intT, n);
    for (intT i=0; i < n; i++) I[i] = i;
    intT m = serialQuickHull(I, P, n);
    cout << "serial-hull-time = " << t.stop() << endl;

    check(P, n, I, m);
    return _seq<intT>(I, m);
  }

  if (sortPoint) {
    auto xLess = [&](point2d p1, point2d p2) {
		   return p1.x()<p2.x();};
    sampleSort(P, n, xLess);
    cout << "sort-time = " << t1.next() << endl;
  }

  intT blk = n/procs;
  intT M[procs+1];

  intT* I = newA(intT, n);

  parallel_for(0, procs,
	       [&](intT i) {
		 intT o = blk*i;
		 if (i != procs-1) {
		   for(intT j=0; j<blk; ++j) I[j+o] = j;
		   M[i] = serialQuickHull(I+o, P+o, blk);
		   //check(P+o, blk, I+o, M[i]);
		 } else {
		   intT blkLast = max(blk, n-blk*(procs-1));
		   for(intT j=0; j<blkLast; ++j) I[j+o] = j;
		   M[i] = serialQuickHull(I+o, P+o, blkLast);
		   //check(P+o, blkLast, I+o, M[i]);
		 }
	       }, 1);
  cout << "subhull-time = " << t1.next() << endl;

  for(intT i=1; i<procs; ++i) {
    M[i] += M[i-1];}
  for(intT i=procs; i>0; --i) {
    M[i] = M[i-1];}
  M[0] = 0;
  cout << "subhull-size = " << M[procs] << endl;

  point2d* PP = newA(point2d, M[procs]);

  parallel_for(0, procs, [&](intT i) {
			   intT o = blk*i;
			   for(intT j=M[i]; j<M[i+1]; ++j) {
			     PP[j] = P[o+I[o+(j-M[i])]];}
			 }, 1);

  intT m2;
  if (M[procs] < 2000) {
    parallel_for(0, M[procs], [&](intT i) {I[i] = i;});
    m2 = serialQuickHull(I, PP, M[procs]);
  } else {
    auto CH = hullInternal(PP, M[procs]);
    m2 = CH.n;
    free(I);
    I = CH.A;
  }
  cout << "merge-time = " << t1.next() << endl;

  cout << "total-hull-time = " << t.stop() << endl;

  //hull is in (PP, I), size is m2
  check(P, n, I, m2, PP);

  free(I);
  free(PP);
  return _seq<intT>();
}
