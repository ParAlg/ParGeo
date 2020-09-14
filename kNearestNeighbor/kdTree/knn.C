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
point<dim>** knn(point<dim>* P, intT n, intT k=1) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;

  static const bool serial = false;
  cout << k << "-nearest neighbor, " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;

  t0.start();
  treeT* tree = new treeT(P, n, !serial, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  pointT** A = newA(pointT*, k*n);
  par_for (intT i=0; i<n; ++i) {
    pointT** NN = tree->kNN(&P[i], k, A+i*k);//query call
  }
  cout << "knn-time = " << t0.next() << endl;
  return A;
}

template point<2>** knn<2>(point<2>*, intT, intT);
template point<3>** knn<3>(point<3>*, intT, intT);
template point<4>** knn<4>(point<4>*, intT, intT);
template point<5>** knn<5>(point<5>*, intT, intT);
template point<6>** knn<6>(point<6>*, intT, intT);
template point<7>** knn<7>(point<7>*, intT, intT);
template point<8>** knn<8>(point<8>*, intT, intT);
template point<9>** knn<9>(point<9>*, intT, intT);
