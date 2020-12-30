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
void perfTest(point<dim>* P, intT n, intT N=-1) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;

  if (N<0) N=n;

  timing t0;

  intT test = 0;
  double times[20];

  auto randIndices = [&](intT n) {
		       intT* I = newA(intT, n);
		       parallel_for(0, n, [&](intT i) {I[i] = i;});
		       randPerm(I, n);
		       return I;
		     };

  {
    intT* A = newA(intT, N);
    intT* B = newA(intT, N);
    t0.start();
    parallel_for(0, N, [&](intT i) {
			 A[i] = B[i];});
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - seq-read: " << times[test++] << " sec" << endl;
  }

  {
    intT* A = newA(intT, N);
    intT* B = newA(intT, N);
    t0.start();
    parallel_for(0, N, [&](intT i) {
			 A[i] = B[i];});
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - seq-write: " << times[test++] << " sec" << endl;
  }

  {
    intT* A = newA(intT, N);
    intT* B = newA(intT, N);
    intT* I = randIndices(N);
    t0.start();
    parallel_for(0, N, [&](intT i) {
			 A[i] = B[I[i]];});
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - random-read: " << times[test++] << " sec" << endl;
  }

  {
    intT* A = newA(intT, N);
    intT* B = newA(intT, N);
    intT* I = randIndices(N);
    t0.start();
    parallel_for(0, N, [&](intT i) {
			 A[I[i]] = B[i];});
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - random-write: " << times[test++] << " sec" << endl;
  }

  {
    t0.start();
    auto I = randIndices(N);
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - rand-perm: " << times[test++] << " sec" << endl;
  }

  {
    intT* A = newA(intT, N);
    intT* I = randIndices(N);
    parallel_for(0, N, [&](intT i) {
			 A[i] = I[i];});
    t0.start();
    sampleSort(A, N, std::less<intT>());
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - sample-sort: " << times[test++] << endl;
  }

  {
    intT* B = newA(intT, N);
    parallel_for(0, N, [&](intT i) {
			 B[i] = utils::hash(i)%2;});
    t0.start();
    intT tmp = sequence::prefixSum(B, 0, N);
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - prefix-sum: " << times[test++] << " sec" << endl;
  }

  {
    t0.start();
    treeT* tree = new treeT(P, n, true, 1);
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - kdtree-construction: " << times[test++] << " sec" << endl;

    intT k = 2;
    pointT** AA = newA(pointT*, k*n);
    t0.start();
    parallel_for (0, n,[&](intT i) {
			 pointT** NN = tree->kNN(&P[i], k, AA+i*k);});
    times[test] = t0.stop();
    cout << ">>> Test " << test << " - knn: " << times[test++] << " sec" << endl;
  }

  cout << ">>>>>>>>" << endl;
  for (intT i=0; i<test; ++i) cout << times[i] << endl;

}

#endif
