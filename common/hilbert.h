#ifndef HILBERT_H
#define HILBERT_H

#include <vector>
#include <iostream>
#include "geometry.h"
#include "pbbs/utils.h"
#include "pbbs/sequence.h"

using namespace std;

template<class T, class cmpT>
inline intT splitItemSerial(T* items, intT n, cmpT leftSide) {
  if (n <= 1) {// todo check
    return 1;
  }
  intT lPt = 0;
  intT rPt = n-1;
  while (lPt < rPt) {
    if (!leftSide(items[lPt])) {
      while (!leftSide(items[rPt]) && lPt < rPt) {
        rPt--;
      }
      if (lPt < rPt) {
        swap(items[lPt], items[rPt]);
        rPt--; }
      else { break;}
    }
    lPt++;
  }
  if (leftSide(items[lPt])) lPt++;
  return lPt;
}

template<class T, class cmpT>
inline intT splitItem(T* A, intT n, cmpT leftSide, T* B, intT* flag) {
  if (n <= 1) {// todo check
    return 1;
  }
  if (n < 2000) return splitItemSerial(A, n, leftSide);
  par_for(intT i=0; i<n; ++i) {
    if (leftSide(A[i])) flag[i]=1;
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

template<class T>
void hilbertSplitMiddleSerial(T* A, intT n, int axe, bool orient, floatT value, intT& middle) {
  if (n <= 0) middle = n;

  auto cmp = [&](T item)
    {
     return orient ? (item[axe] > value) : (item[axe] <= value);
    };
  middle = splitItemSerial(A, n, cmp);
}

template<class T>
void hilbertSplitMiddle(T* A, intT s, intT e, int axe, bool orient, floatT value, intT& middle, T* B, intT* flag) {
  if (e-s <= 0) middle = e-s;

  auto cmp = [&](T item)
    {
     return orient ? (item[axe] > value) : (item[axe] <= value);
    };
  middle = splitItem(A+s, e-s, cmp, B, flag);
  middle += s;//offset
}

template<int dim, class T>
void hilbertMiddleHelperSerial(T* A, intT n, bool* startt, intT direction, point<dim> pMin, point<dim> pMax, intT numQuadrant) {
  if (n <= 1) return;

  typedef point<dim> pointT;

  pointT med = pMin.average(pMax);
  pointT cmin = pMin;
  pointT cmax = med;

  bool start[dim]; for(int i=0; i<dim; ++i) start[i]=startt[i];
  int places[numQuadrant +1];
  int dir[numQuadrant +1];
  places[0] = 0;
  places[numQuadrant] = n;

  int lastDir = (direction + dim) % dim;
  int curDir = direction;
  int levelStep = numQuadrant;
  do{
    int halfStep = levelStep/2;
    int left = 0;
    int middle = halfStep;
    int right = levelStep;
    bool orient = start[curDir];

    do{
      dir[middle] = curDir;
      hilbertSplitMiddleSerial(A+places[left], places[right]-places[left], curDir, orient, med[curDir], places[middle]);
      places[middle] += places[left];//offset
      left = right;
      right += levelStep;
      middle += levelStep;
      orient = !orient;
    } while (left < numQuadrant);

    levelStep = halfStep;
    curDir = (curDir+1) % dim;
  } while (curDir != lastDir);

  lastDir = (direction + dim -1) % dim;

  if (places[1] != n)
    hilbertMiddleHelperSerial(A+places[0], places[1]-places[0], start, lastDir, cmin, cmax, numQuadrant);

  cmin[lastDir] = med[lastDir];
  cmax[lastDir] = pMax[lastDir];

  for(int i=1; i<numQuadrant-1; i+=2){

    if (places[i]!=0 || places[i+1]!=n)
      hilbertMiddleHelperSerial(A+places[i], places[i+1]-places[i], start, dir[i+1], cmin, cmax, numQuadrant);

    cmax[dir[i+1]] = (cmin[dir[i+1]] == pMin[dir[i+1]])
      ? pMax[dir[i+1]] : pMin[dir[i+1]];
    cmin[dir[i+1]] = med[dir[i+1]];

    if (places[i+1]!=0 || places[i+2]!=n)
      hilbertMiddleHelperSerial(A+places[i+1], places[i+2]-places[i+1], start, dir[i+1], cmin, cmax, numQuadrant);

    cmin[dir[i+1]] = cmax[dir[i+1]];
    cmax[dir[i+1]] = med[dir[i+1]];
    cmax[lastDir] = (cmax[lastDir]==pMax[lastDir])
      ? pMin[lastDir] : pMax[lastDir];

    start[dir[i+1]] = !start[dir[i+1]];
    start[lastDir] = !start[lastDir];
  }

  if (places[numQuadrant-1]!=0)
    hilbertMiddleHelperSerial(A+places[numQuadrant-1], places[numQuadrant]-places[numQuadrant-1], start, lastDir, cmin, cmax, numQuadrant);
}

template<int dim, class T>
void hilbertMiddleHelper(T* A, intT n, bool* startt, intT direction, point<dim> pMin, point<dim> pMax, intT numQuadrant, T* B, intT* flag) {
  if (n <= 1) return;
  if (n < 2000) return hilbertMiddleHelperSerial(A, n, startt, direction, pMin, pMax, numQuadrant);

  static const bool verbose = false;

  typedef point<dim> pointT;

  pointT med = pMin.average(pMax);
  pointT cmin = pMin;
  pointT cmax = med;

  bool start[dim]; for(int i=0; i<dim; ++i) start[i]=startt[i];
  int places[numQuadrant +1];
  int dir[numQuadrant +1];
  places[0] = 0;
  places[numQuadrant] = n;

  int lastDir = (direction + dim) % dim;
  int curDir = direction;
  int levelStep =numQuadrant;
  do{
    int halfStep = levelStep/2;
    int left = 0;
    int middle = halfStep;
    int right = levelStep;
    bool orient = start[curDir];

    do{
      dir[middle] = curDir;

      if (verbose) {
        cout << "splitting A[" << places[left] << ":" << places[right] << "]" << endl;
        cout << " split, n = " << places[right]-places[left] << ", orient = " << orient << ", axe = " << curDir << ", value = " << med[curDir] << endl;
      }

      cilk_spawn hilbertSplitMiddle(A, places[left], places[right],
                                    curDir, orient, med[curDir], places[middle],
                                    B+places[left], flag+places[left]);

      left = right;
      right += levelStep;
      middle += levelStep;
      orient = !orient;
    } while (left < numQuadrant);
    cilk_sync;

    levelStep = halfStep;
    curDir = (curDir+1) % dim;
  } while (curDir != lastDir);

  lastDir = (direction + dim -1) % dim;

  if (places[1] != n) {
    cilk_spawn hilbertMiddleHelper(A+places[0], places[1]-places[0], start, lastDir, cmin, cmax, numQuadrant, B+places[0], flag+places[0]);
  }

  cmin[lastDir] = med[lastDir];
  cmax[lastDir] = pMax[lastDir];

  for(int i=1; i<numQuadrant-1; i+=2){

    if (places[i]!=0 || places[i+1]!=n) {
      cilk_spawn hilbertMiddleHelper(A+places[i], places[i+1]-places[i], start, dir[i+1], cmin, cmax, numQuadrant, B+places[i], flag+places[i]);
    }

    cmax[dir[i+1]] = (cmin[dir[i+1]] == pMin[dir[i+1]])
      ? pMax[dir[i+1]] : pMin[dir[i+1]];
    cmin[dir[i+1]] = med[dir[i+1]];

    if (places[i+1]!=0 || places[i+2]!=n) {
      cilk_spawn hilbertMiddleHelper(A+places[i+1], places[i+2]-places[i+1], start, dir[i+1], cmin, cmax, numQuadrant, B+places[i+1], flag+places[i+1]);
    }

    cmin[dir[i+1]] = cmax[dir[i+1]];
    cmax[dir[i+1]] = med[dir[i+1]];
    cmax[lastDir] = (cmax[lastDir]==pMax[lastDir])
      ? pMin[lastDir] : pMax[lastDir];

    start[dir[i+1]] = !start[dir[i+1]];
    start[lastDir] = !start[lastDir];
  }

  if (places[numQuadrant-1]!=0) {
    cilk_spawn hilbertMiddleHelper(A+places[numQuadrant-1], places[numQuadrant]-places[numQuadrant-1], start, lastDir, cmin, cmax, numQuadrant, B+places[numQuadrant-1], flag+places[numQuadrant-1]);
  }
  cilk_sync;
}

template<int dim, class T>
void hilbertMiddle(T* A, intT n) {
  typedef point<dim> pointT;

  auto box = boundingBoxParallel<dim, T>(A, n);
  auto pMin = box.first;
  auto pMax = box.second;

  intT numQuadrant = 1;
  bool start[dim];
  for (int i=0; i<dim; ++i) {
    start[i]=false;
    numQuadrant *= 2;
  }

  if (n < 2000) {
    hilbertMiddleHelperSerial<dim, T>(A, n, start, 0, pMin, pMax, numQuadrant);
  } else {
    auto B = newA(T, n);
    auto flag = newA(intT, n);
    hilbertMiddleHelper<dim, T>(A, n, start, 0, pMin, pMax, numQuadrant, B, flag);
    free(B);
    free(flag);
  }
}

#endif
