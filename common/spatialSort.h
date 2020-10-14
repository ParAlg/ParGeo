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

#ifndef SPATIAL_SORT_H
#define SPATIAL_SORT_H

#include "geometry.h"
#include "pbbs/sequence.h"
#include "pbbs/utils.h"

template<int dim, class T>
pair<point<dim>, point<dim>> boundingBoxSerial(T* A, intT n) {
  typedef point<dim> pointT;
  auto pMin = pointT(A[0].coordinate());
  auto pMax = pointT(A[0].coordinate());
  for(intT i=0; i<n; ++i) {
    pMin.minCoords(A[i].coordinate());
    pMax.maxCoords(A[i].coordinate());
  }
  return make_pair(pMin, pMax);
}

template<int dim, class T>
pair<point<dim>, point<dim>> boundingBoxParallel(T* A, intT n) {
  typedef point<dim> pointT;
  intT P = getWorkers()*8;
  intT blockSize = (n+P-1)/P;
  pointT localMin[P];
  pointT localMax[P];
  for (intT i=0; i<P; ++i) {
    localMin[i] = pointT(A[0].coordinate());
    localMax[i] = pointT(A[0].coordinate());}
  par_for(intT p=0; p<P; ++p) {
    intT s = p*blockSize;
    intT e = min((intT)(p+1)*blockSize,n);
    for (intT j=s; j<e; ++j) {
      localMin[p].minCoords(A[j].coordinate());
      localMax[p].maxCoords(A[j].coordinate());}
  }
  auto pMin = pointT(A[0].coordinate());
  auto pMax = pointT(A[0].coordinate());
  for(intT p=0; p<P; ++p) {
    pMin.minCoords(localMin[p].x);
    pMax.maxCoords(localMax[p].x);}
  return make_pair(pMin, pMax);
}

template<int dim, class T>
inline intT splitItemSerial(T* A, intT n, floatT xM, intT k) {
  if (n < 2) {
    cout << "error, spatial sort splitting singleton, abort" << endl;abort();}
  intT lPt = 0;
  intT rPt = n-1;
  while (lPt < rPt) {
    if (A[lPt].coordinate(k)>=xM) {
      while (A[rPt].coordinate(k)>=xM && lPt < rPt) {
        rPt--;
      }
      if (lPt < rPt) {
        swap(A[lPt], A[rPt]);
        rPt--; }
      else { break;}
    }
    lPt++;
  }
  if (A[lPt].coordinate(k) < xM) lPt++;
  return lPt;
}

template<int dim, class T>
inline intT splitItemParallel(T* A, intT n, floatT xM, intT k, T* B, intT* flag) {
  if (n < 2) {
    cout << "error, spatial sort splitting singleton, abort" << endl;abort();}
  par_for(intT i=0; i<n; ++i) {
    if (A[i].coordinate(k)<xM) flag[i]=1;
    else flag[i] = 0;
  }
  intT leftSize = sequence::prefixSum(flag,0,n);
  par_for(intT i=0; i<n-1; ++i) {
    if (flag[i] != flag[i+1]) B[flag[i]] = A[i];
    if (i-flag[i] != i+1-flag[i+1]) B[leftSize+i-flag[i]] = A[i];
  }
  if (flag[n-1] != leftSize) B[flag[n-1]] = A[n-1];
  if (n-1-flag[n-1] != n-leftSize) B[leftSize+n-1-flag[n-1]] = A[n-1];
  par_for(intT i=0; i<n; ++i) {
    A[i] = B[i];
  }
  return leftSize;
}

// T - a spatial data type of k-dim, supports coordinate(i) to get i-th dim
// cmpT - a spatial comparator that can take in dimension
template<int dim, class T>
void spatialSortSerial(T* A, intT n, intT thresh=16) {
  if (n <= thresh) return;

  auto bb = boundingBoxSerial<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  floatT xM = -1;
  intT k;
  for (int kk=0; kk<dim; ++kk) {
    if (pMax[kk]-pMin[kk]>xM) {
      xM = pMax[kk]-pMin[kk];
      k = kk;}}
  xM = (pMax[k]+pMin[k])/2;
  intT median = splitItemSerial<dim, T>(A, n, xM, k);
  spatialSortSerial<dim, T>(A, median, thresh);
  spatialSortSerial<dim, T>(A+median, n-median, thresh);
}

template<int dim, class T>
void spatialSort(T* A, intT n, intT thresh=16, T* B=NULL, intT* flag=NULL) {
  if (n <= thresh) return;
  if (n < 2000) return spatialSortSerial<dim, T>(A, n);

  bool freeB = false;
  if (!B) {
    B = newA(T, n);
    freeB = true;}

  bool freeFlag = false;
  if (!flag) {
    flag = newA(intT, n);
    freeFlag = true;}

  auto bb = boundingBoxParallel<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  floatT xM = -1;
  intT k;
  for (int kk=0; kk<dim; ++kk) {
    if (pMax[kk]-pMin[kk]>xM) {
      xM = pMax[kk]-pMin[kk];
      k = kk;}}
  xM = (pMax[k]+pMin[k])/2;
  intT median = splitItemParallel<dim, T>(A, n, xM, k, B, flag);

  cilk_spawn spatialSort<dim, T>(A, median, thresh, B, flag);
  spatialSort<dim, T>(A+median, n-median, thresh, B+median, flag+median);
  cilk_sync;

  if(freeB) free(B);
  if(freeFlag) free(flag);
}

#endif
