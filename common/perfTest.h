#ifndef PERF_TEST_H
#define PERF_TEST_H

#include "geometry.h"
#include "kdTree.h"
#include "kNearestNeighbors.h"
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "pbbs/sampleSort.h"

template <int dim>
void perfTest(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;

  timing t0;

  intT* A = newA(intT, n);
  intT* B = newA(intT, n);
  double x1, x2, x3, x4, x5;

  t0.start();
  parallel_for(0, n, [&](intT i) {
      A[i] = utils::hash(i);
      B[i] = utils::hash(i)%2;
    });
  x1 = t0.next();
  cout << ">>> Test 1 - par-for IO: " << x1 << " sec" << endl;

  intT tmp = sequence::prefixSum(B, 0, n);
  x2 = t0.next();
  cout << ">>> Test 2 - prefixSum: " << x2 << " sec" << endl;

  sampleSort(A, n, std::less<intT>());
  x3 = t0.next();
  cout << ">>> Test 3 - sampleSort: " << x3 << endl;

  treeT* tree = new treeT(P, n, true, 1);
  x4 = t0.next();
  cout << ">>> Test 4 - kdTree construction: " << x4 << " sec" << endl;

  intT k = 2;
  pointT** AA = newA(pointT*, k*n);
  t0.next();
  parallel_for (0, n,[&](intT i) {
		       pointT** NN = tree->kNN(&P[i], k, AA+i*k);});
  x5 = t0.next();
  cout << ">>> Test 5 - data parallel knn: " << x5 << " sec" << endl;

  cout << ">>>>>>>>" << endl;
  cout << x1 << endl;
  cout << x2 << endl;
  cout << x3 << endl;
  cout << x4 << endl;
  cout << x5 << endl;
}

#endif
