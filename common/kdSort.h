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
#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"
#include <iterator>
#include <algorithm>
#include <queue>

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

template<int dim>
intT findWidest(point<dim> pMin, point<dim> pMax) {
  floatT xM = -1;
  intT k;
  for (int kk=0; kk<dim; ++kk) {
    if (pMax[kk]-pMin[kk]>xM) {
      xM = pMax[kk]-pMin[kk];
      k = kk;}}
  return k;
}

template<int dim, class T>
intT findMiddleSerial(T* A, intT n, point<dim> pMin, point<dim> pMax) {
  auto k = findWidest(pMin, pMax);
  floatT xM = (pMax[k]+pMin[k])/2;
  intT median = splitItemSerial<dim, T>(A, n, xM, k);
  return median;
}

template<int dim, class T>
intT findMiddleParallel(T* A, intT n, point<dim> pMin, point<dim> pMax, T* B, intT* flag) {
  auto k = findWidest(pMin, pMax);
  floatT xM = (pMax[k]+pMin[k])/2;
  intT median = splitItemParallel<dim, T>(A, n, xM, k, B, flag);
  return median;
}

// T - a spatial data type of k-dim, supports coordinate(i) to get i-th dim
// cmpT - a spatial comparator that can take in dimension
template<int dim, class T>
void kdSortMiddleSerial(T* A, intT n, intT thresh=16) {
  if (n <= thresh) return;

  auto bb = boundingBoxSerial<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  intT median = findMiddleSerial(A, n, pMin, pMax);
  if (median == 0 || median == n) {median = ceil(n/2.0);}

  kdSortMiddleSerial<dim, T>(A, median, thresh);
  kdSortMiddleSerial<dim, T>(A+median, n-median, thresh);
}

template<int dim, class T>
void kdSortMiddle(T* A, intT n, intT thresh=16, T* B=NULL, intT* flag=NULL) {
  if (n <= thresh) return;
  if (n < 2000) return kdSortMiddleSerial<dim, T>(A, n);

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

  intT median = findMiddleParallel(A, n, pMin, pMax, B, flag);
  if (median == 0 || median == n) {median = ceil(n/2.0);}

  cilk_spawn kdSortMiddle<dim, T>(A, median, thresh, B, flag);
  kdSortMiddle<dim, T>(A+median, n-median, thresh, B+median, flag+median);
  cilk_sync;

  if(freeB) free(B);
  if(freeFlag) free(flag);
}

template<int dim, class T>
void kdSortMedianSerial(T* A, intT n, intT thresh=16) {
  if (n <= thresh) return;
  auto bb = boundingBoxSerial<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  auto k = findWidest(pMin, pMax);

  auto splitK = [&](T a, T b) {
                  return a[k] < b[k];
                };

  intT median = ceil(n/2.0);
  std::nth_element(A, A+median, A+n, splitK);

  kdSortMedianSerial<dim, T>(A, median, thresh);
  kdSortMedianSerial<dim, T>(A+median, n-median, thresh);
}

template<int dim, class T>
void kdSortMedian(T* A, intT n, intT thresh=16) {
  if (n <= thresh) return;
  if (n <= 2000) return kdSortMedianSerial<dim, T>(A, n, thresh);

  auto bb = boundingBoxParallel<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  auto k = findWidest(pMin, pMax);

  auto splitK = [&](T a, T b) {
                  return a[k] < b[k];
                };

  intT median = ceil(n/2.0);
  std::nth_element(A, A+median, A+n, splitK);//todo parallel?
  //sampleSort(A, n, splitK);

  cilk_spawn kdSortMedian<dim, T>(A, median, thresh);
  kdSortMedian<dim, T>(A+median, n-median, thresh);
  cilk_sync;
}

template<class T>
struct bfsNode {
  T item;
  bfsNode* left;
  bfsNode* right;
  bool visited;
};

template<int dim, class T>
void kdSortBFSHelperSerial(T* A, intT n, bfsNode<T>* space, bfsNode<T>* start, intT thresh=16) {
  if (n <= thresh) {
    //mark as leaf
    space[0].item = A[0];
    space->left = NULL;
    space->right = space->left + n;
    for (intT i=1; i<n; ++i) {
      // if (!space[i].item.isEmpty()) {
      //   cout << "error, overwrite " << space-start << endl;
      //   abort();
      // }
      space[i].item = A[i];
      space[i].left = NULL;
      space[i].right = NULL;
    }
    return;
  }
  // if (!space->item.isEmpty()) {
  //   cout << "error, overwrite " << space-start << endl;
  //   abort();
  // }
  auto bb = boundingBoxSerial<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  auto k = findWidest(pMin, pMax);

  auto splitK = [&](T a, T b) {
                  return a[k] < b[k];
                };

  intT median = ceil(n/2.0);
  //cout << "median = " << median << " out of " << n << ", write to " << space-start << endl;
  std::nth_element(A, A+median, A+n, splitK);
  if (median == 0 || median == n) {median = n/2;}
  space->item = A[median-1];
  space->left = space+1;
  space->right = space+median;

  kdSortBFSHelperSerial<dim, T>(A, median-1, space->left, start, thresh);
  kdSortBFSHelperSerial<dim, T>(A+median, n-median, space->right, start, thresh);
}

template<int dim, class T>
void kdSortBFSHelperParallel(T* A, intT n, bfsNode<T>* space, bfsNode<T>* start, intT thresh=16) {
  if (n <= thresh) {
    //mark as leaf
    space[0].item = A[0];
    space->left = NULL;
    space->right = space->left + n;
    for (intT i=1; i<n; ++i) {
      // if (!space[i].item.isEmpty()) {
      //   cout << "error, overwrite " << space-start << endl;
      //   abort();
      // }
      space[i].item = A[i];
      space[i].left = NULL;
      space[i].right = NULL;
    }
    return;
  }
  // if (!space->item.isEmpty()) {
  //   cout << "error, overwrite " << space-start << endl;
  //   abort();
  // }

  if (n < 2000) return kdSortBFSHelperSerial<dim, T>(A, n, space, start, thresh);

  auto bb = boundingBoxParallel<dim, T>(A, n);
  auto pMin = bb.first; auto pMax = bb.second;
  auto k = findWidest(pMin, pMax);

  auto splitK = [&](T a, T b) {
                  return a[k] < b[k];
                };

  intT median = ceil(n/2.0);
  //cout << "median = " << median << " out of " << n << ", write to " << space-start << endl;
  std::nth_element(A, A+median, A+n, splitK);
  if (median == 0 || median == n) {median = n/2;}
  space->item = A[median-1];
  space->left = space+1;
  space->right = space+median;

  cilk_spawn kdSortBFSHelperParallel<dim, T>(A, median-1, space->left, start, thresh);
  kdSortBFSHelperParallel<dim, T>(A+median, n-median, space->right, start, thresh);
  cilk_sync;
}

template<class T>
void kdBFS(bfsNode<T>* root, T* write) {
  auto isLeaf = [&](bfsNode<T>* node) {
                  if (!node->left && node->right) return true;
                  else return false;
                };
  auto leafSize = [&](bfsNode<T>* node) {
                    return node->right - node->left;
                  };

  auto Q = queue<bfsNode<T>*>();
  Q.push(root);
  root->visited = true;

  intT ii=0;
  while (Q.size() > 0) {
    auto visit = Q.front();
    Q.pop();
    write[ii++] = visit->item;

    if (isLeaf(visit)) {
      for (intT i=1; i<visit->right-visit->left; ++i) {
        write[ii++] = visit[i].item;
      }
    } else {
      if(!visit->left->visited) {
        visit->left->visited = true;
        Q.push(visit->left);}
      if(!visit->right->visited) {
        visit->right->visited = true;
        Q.push(visit->right);}
    }
  }
}

template<int dim, class T>
void kdSortBFSSerial(T* A, intT n, intT thresh=16) {
  auto space = newA(bfsNode<T>, n);
  for (int i=0; i<n; ++i) {
    space[i].visited = false;
    //space[i].item = T();//empty
  }
  kdSortBFSHelperSerial<dim>(A, n, space, space, thresh);
  kdBFS<T>(space, A);
  free(space);
}

template<int dim, class T>
void kdSortBFS(T* A, intT n, intT thresh=16) {
  auto space = newA(bfsNode<T>, n);
  for (int i=0; i<n; ++i) {
    space[i].visited = false;
    //space[i].item = T();//empty
  }
  kdSortBFSHelperParallel<dim>(A, n, space, space, thresh);
  timing t; t.start();
  kdBFS<T>(space, A);
  cout << "bfs-time(seq) = " << t.stop() << endl;
  free(space);
}

#endif
