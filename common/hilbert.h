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
void hilbertSplitMiddle(T* A, intT n, int axe, bool orient, floatT value, intT& middle, T* B, intT* flag) {
  if (n <= 0) middle = n;

  auto cmp = [&](T item)
    {
     return orient ? (item[axe] > value) : (item[axe] <= value);
    };
  middle = splitItem(A, n, cmp, B, flag);
}

template<int dim, class T>
void hilbertMiddleHelperSerial(T* A, intT n, vector<bool> start, intT direction, point<dim> mini, point<dim> maxi, intT two_to_dim) {
  if (n <= 1) return;

  typedef point<dim> pointT;

  pointT med = mini.average(maxi);
  pointT cmin = mini;
  pointT cmax = med;

  std::vector<int> places(two_to_dim +1);
  std::vector<int> dir(two_to_dim +1);
  places[0] = 0;
  places[two_to_dim] = n;

  int last_dir = (direction + dim) % dim;
  int current_dir = direction;
  int current_level_step =two_to_dim;
  do{
    int half_step = current_level_step/2;
    int left = 0;
    int middle = half_step;
    int right = current_level_step;
    bool orient = start[current_dir];

    do{
      dir[middle] = current_dir;

      hilbertSplitMiddleSerial(A+places[left], places[right]-places[left], current_dir, orient, med[current_dir], places[middle]);
      places[middle] += places[left];//offset

      left = right;
      right += current_level_step;
      middle += current_level_step;
      orient = !orient;
    } while (left < two_to_dim);

    current_level_step = half_step;
    current_dir = (current_dir+1) % dim;
  } while (current_dir != last_dir);

  last_dir = (direction + dim -1) % dim;

  if (places[1] != n)
    hilbertMiddleHelperSerial(A+places[0], places[1]-places[0], start, last_dir, cmin, cmax, two_to_dim);

  cmin[last_dir] = med[last_dir];
  cmax[last_dir] = maxi[last_dir];

  for(int i=1; i<two_to_dim-1; i+=2){

    if (places[i]!=0 || places[i+1]!=n)
      hilbertMiddleHelperSerial(A+places[i], places[i+1]-places[i], start, dir[i+1], cmin, cmax, two_to_dim);

    cmax[dir[i+1]] = (cmin[dir[i+1]] == mini[dir[i+1]])
      ? maxi[dir[i+1]] : mini[dir[i+1]];
    cmin[dir[i+1]] = med[dir[i+1]];

    if (places[i+1]!=0 || places[i+2]!=n)
      hilbertMiddleHelperSerial(A+places[i+1], places[i+2]-places[i+1], start, dir[i+1], cmin, cmax, two_to_dim);

    cmin[dir[i+1]] = cmax[dir[i+1]];
    cmax[dir[i+1]] = med[dir[i+1]];
    cmax[last_dir] = (cmax[last_dir]==maxi[last_dir])
      ? mini[last_dir] : maxi[last_dir];

    start[dir[i+1]] = !start[dir[i+1]];
    start[last_dir] = !start[last_dir];
  }

  if (places[two_to_dim-1]!=0)
    hilbertMiddleHelperSerial(A+places[two_to_dim-1], places[two_to_dim]-places[two_to_dim-1], start, last_dir, cmin, cmax, two_to_dim);
}

template<int dim, class T>
void hilbertMiddleHelper(T* A, intT n, vector<bool> start, intT direction, point<dim> mini, point<dim> maxi, intT two_to_dim, T* B, intT* flag) {
  if (n < 2000) return hilbertMiddleHelperSerial(A, n, start, direction, mini, maxi, two_to_dim);

  static const bool verbose = false;

  if (n <= 1) return;

  typedef point<dim> pointT;

  pointT med = mini.average(maxi);
  pointT cmin = mini;
  pointT cmax = med;

  std::vector<int> places(two_to_dim +1);
  std::vector<int> dir(two_to_dim +1);
  places[0] = 0;
  places[two_to_dim] = n;

  int last_dir = (direction + dim) % dim;
  int current_dir = direction;
  int current_level_step =two_to_dim;
  do{
    int half_step = current_level_step/2;
    int left = 0;
    int middle = half_step;
    int right = current_level_step;
    bool orient = start[current_dir];

    do{
      dir[middle] = current_dir;

      if (verbose) {
        cout << "splitting A[" << places[left] << ":" << places[right] << "]" << endl;
        cout << " split, n = " << places[right]-places[left] << ", orient = " << orient << ", axe = " << current_dir << ", value = " << med[current_dir] << endl;
        cout << " left = " << left << ", middle = " << middle << ", right = " << right << endl;
      }

      hilbertSplitMiddle(A+places[left], places[right]-places[left], current_dir, orient, med[current_dir], places[middle], B, flag);
      places[middle] += places[left];//offset

      if(verbose) {
        cout << "leftSize = " << places[middle]-places[left] << endl;
        cout << "rightSize = " << places[right]-places[middle] << endl;
        //for (intT i=places[left]; i<places[middle]; ++i) cout << A[i] << " ";cout << endl << endl;
        //for (intT i=places[middle]; i<places[right]; ++i) cout << A[i] << " ";cout << endl << endl;
        cout << "--------------------" << endl;
      }

      left = right;
      right += current_level_step;
      middle += current_level_step;
      orient = !orient;
    } while (left < two_to_dim);

    current_level_step = half_step;
    current_dir = (current_dir+1) % dim;
  } while (current_dir != last_dir);

  last_dir = (direction + dim -1) % dim;

  if (places[1] != n) {
    cilk_spawn hilbertMiddleHelper(A+places[0], places[1]-places[0], start, last_dir, cmin, cmax, two_to_dim, B+places[0], flag+places[0]);
  }

  cmin[last_dir] = med[last_dir];
  cmax[last_dir] = maxi[last_dir];

  for(int i=1; i<two_to_dim-1; i+=2){

    if (places[i]!=0 || places[i+1]!=n) {
      cilk_spawn hilbertMiddleHelper(A+places[i], places[i+1]-places[i], start, dir[i+1], cmin, cmax, two_to_dim, B+places[i], flag+places[i]);
    }

    cmax[dir[i+1]] = (cmin[dir[i+1]] == mini[dir[i+1]])
      ? maxi[dir[i+1]] : mini[dir[i+1]];
    cmin[dir[i+1]] = med[dir[i+1]];

    if (places[i+1]!=0 || places[i+2]!=n) {
      cilk_spawn hilbertMiddleHelper(A+places[i+1], places[i+2]-places[i+1], start, dir[i+1], cmin, cmax, two_to_dim, B+places[i+1], flag+places[i+1]);
    }

    cmin[dir[i+1]] = cmax[dir[i+1]];
    cmax[dir[i+1]] = med[dir[i+1]];
    cmax[last_dir] = (cmax[last_dir]==maxi[last_dir])
      ? mini[last_dir] : maxi[last_dir];

    if(verbose) {
      cout << "start[" << dir[i+1] << "] = " << !start[dir[i+1]] << endl;
      cout << "start[" << last_dir << "] = " << !start[last_dir] << endl;}
    start[dir[i+1]] = !start[dir[i+1]];
    start[last_dir] = !start[last_dir];
  }

  if (places[two_to_dim-1]!=0) {
    if(verbose) cout << "last step" << endl;
    cilk_spawn hilbertMiddleHelper(A+places[two_to_dim-1], places[two_to_dim]-places[two_to_dim-1], start, last_dir, cmin, cmax, two_to_dim, B+places[two_to_dim-1], flag+places[two_to_dim-1]);
  }
  cilk_sync;
}

template<int dim, class T>
  void hilbertMiddle(T* A, intT n) {
  typedef point<dim> pointT;

  auto box = boundingBoxParallel<dim, T>(A, n);
  auto mini = box.first;
  auto maxi = box.second;

  intT two_to_dim = 1;
  vector<bool> start(dim);
  for (int i=0; i<dim; ++i) {
    start[i]=false;
    two_to_dim *= 2;
  }

  if (n < 2000) {
    hilbertMiddleHelperSerial<dim, T>(A, n, start, 0, mini, maxi, two_to_dim);
  } else {
    auto B = newA(T, n);
    auto flag = newA(intT, n);
    hilbertMiddleHelper<dim, T>(A, n, start, 0, mini, maxi, two_to_dim, B, flag);
    free(B);
    free(flag);
  }
}

#endif
