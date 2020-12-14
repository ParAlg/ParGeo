#ifndef PERF_TEST_H
#define PERF_TEST_H

#include "geometry.h"
#include "kdTree.h"
#include "kNearestNeighbors.h"
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "pbbs/sampleSort.h"
#include "pbbs/randPerm.h"

template <int dim>
void perfTest(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;

  timing t0;

  intT test = 0;
  intT* I = newA(intT, n);
  intT* A = newA(intT, n);
  intT* B = newA(intT, n);
  double times[20];

  t0.start();

  parallel_for(0, n, [&](intT i) {I[i] = i;});
  randPerm(I, n);
  times[test] = t0.next();
  cout << ">>> Test " << test << " - rand-perm: " << times[test++] << " sec" << endl;

  parallel_for(0, n, [&](intT i) {
		       A[I[i]] = i;
		     });
  times[test] = t0.next();
  cout << ">>> Test " << test << " - random-write: " << times[test++] << " sec" << endl;

  parallel_for(0, n, [&](intT i) {
		       A[i] = i;
		     });
  times[test] = t0.next();
  cout << ">>> Test " << test << " - seq-write: " << times[test++] << " sec" << endl;

  parallel_for(0, n, [&](intT i) {
		       auto tmp = A[I[i]];
		     });
  times[test] = t0.next();
  cout << ">>> Test " << test << " - random-read: " << times[test++] << " sec" << endl;

  parallel_for(0, n, [&](intT i) {
		       auto tmp = A[i];
		     });
  times[test] = t0.next();

  cout << ">>> Test " << test << " - seq-read: " << times[test++] << " sec" << endl;

  parallel_for(0, n, [&](intT i) {
		       B[i] = utils::hash(i)%2;
		     });
  t0.next();
  intT tmp = sequence::prefixSum(B, 0, n);
  times[test] = t0.next();
  cout << ">>> Test " << test << " - prefix-sum: " << times[test++] << " sec" << endl;

  sampleSort(A, n, std::less<intT>());
  times[test] = t0.next();
  cout << ">>> Test " << test << " - sample-sort: " << times[test++] << endl;

  treeT* tree = new treeT(P, n, true, 1);
  times[test] = t0.next();
  cout << ">>> Test " << test << " - kdtree-construction: " << times[test++] << " sec" << endl;

  intT k = 2;
  pointT** AA = newA(pointT*, k*n);
  times[test] = t0.next();
  parallel_for (0, n,[&](intT i) {
		       pointT** NN = tree->kNN(&P[i], k, AA+i*k);});
  times[test] = t0.next();
  cout << ">>> Test " << test << " - knn: " << times[test++] << " sec" << endl;

  cout << ">>>>>>>>" << endl;
  for (intT i=0; i<test; ++i) cout << times[i] << endl;

  free(I); free(A); free(B);
}

#endif
