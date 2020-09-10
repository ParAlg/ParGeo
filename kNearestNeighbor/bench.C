#include "kdTree.h"
#include "kdNode.h"
#include "kBuffer.h"
#include "kNearestNeighbors.h"
#include "pbbs/gettime.h"
using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;

  static const bool serial = false;
  cout << "test kd tree, " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;

  t0.start();
  treeT* tree = new treeT(P, n, !serial, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  //query the k-nn of q points
  intT k = 10;
  intT q = 10;
  for (intT i=0; i<q; ++i) {
    pointT** A = tree->kNN(&P[i], k);//query call

    cout << "(i=" << i << ", " << P[i] << ") " << endl;
    for (intT j=0; j<k; ++j) {
      auto nbr = *A[j];
      cout << " " << nbr << ": " << nbr.pointDist(P[i]) << endl;;
    }
    cout << endl;
    free(A);
  }
  cout << "knn-time = " << t0.next() << endl;
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
