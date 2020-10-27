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

#include "pbbs/gettime.h"
#include "pbbs/sampleSort.h"
#include "geometry.h"
#include "kdSort.h"
#include "brio.h"
#include "hilbert.h"
using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  cout << "test spatialSort " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  pointT* Q = newA(pointT, n);
  par_for (int i=0; i<n; ++i) Q[i] = P[i];

  timing t0;
  t0.start();

  //for (intT i=0; i<n; ++i) cout << P[i] << " ";cout << endl << endl;
  hilbertMiddle<dim>(P, n);
  //brioSort<dim>(P, n, n*0.2);
  //kdSortMiddle<dim, pointT>(P, n);
  //kdSortMedian<dim, pointT>(P, n);
  //kdSortBFS<dim, pointT>(P, n);
  //for (intT i=0; i<n; ++i) cout << P[i] << endl;cout << endl << endl;
  //for (intT i=0; i<20; ++i) cout << P[i] << endl;cout << endl << endl;
  //for (intT i=n-20; i<n; ++i) cout << P[i] << endl;cout << endl << endl;
  cout << "spatial-sort-time = " << t0.stop() << endl;

  auto cmp = [&](pointT a, pointT b)
             {
               if (a[0] < b[0]) {
                 return true;
               } else {
                 return false;
               }
             };
  sampleSort(P, n, cmp);
  sampleSort(Q, n, cmp);
  for (int i=0; i<n; ++i) {
    if(Q[i] != P[i]) {
      cout << "something wrong" << endl;
      abort();
    }
  }
  free(Q);
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
