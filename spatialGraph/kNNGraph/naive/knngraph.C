#include "kNearestNeighbor/kdTree/knn.C"
#include "pbbs/gettime.h"
#include "pbbs/sampleSort.h"

#include <cstdint>
using namespace std;

#define DEBUG

// A: chunks of k NNs
template<int dim>
tuple<intT*, intT*> directedKnnGraph(point<dim>* P, point<dim>** A, intT n, intT k, bool sortit){
    typedef point<dim> pointT;
    intT *offsets = (intT*) malloc((n)*sizeof(intT));
    intT *edges = (intT*) malloc((n*k)*sizeof(intT));
    auto base_addr = reinterpret_cast<uintptr_t>(P);
    auto point_size = sizeof(pointT);
    parallel_for (0, n, [&](intT i) {offsets[i] = k;});
    parallel_for (0, n, [&](intT i) {
        parallel_for (0, k, [&](intT j) {
            edges[i*k + j] = (reinterpret_cast<uintptr_t>(A[i*k+j]) - base_addr) / point_size;
        });
    });

    if(sortit){
    auto cmp = [&](intT a, intT b)
        {if (a!=b) return a < b;
         return true;
        };
    parallel_for (0, n, [&](intT i) {
        sampleSort(&edges[i*k], k, cmp);
    });  
    }

// #ifdef DEBUG
//     // for(intT i = 0; i < n; ++i){
//     //     cout << offsets[i] << endl;
//     // }
//     // cout << "======" << endl;
//     for(intT i = 0; i < n; ++i){
//         cout  << ":::" <<  i << ":::" << endl;
//         for(intT j = 0; j < k; ++j){
//             cout << edges[i*k + j] << endl;
//         }
//     }
//     cout << "======" << endl;
// #endif

    return make_tuple(offsets, edges);
}

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
tuple<intT*, intT*> knnGraph(point<dim>* P, intT n, intT k, bool directed, bool sortit) {
  typedef point<dim> pointT;

  cout << k << "-nearest neighbor graph construction, " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;
  pointT** A = knn(P, n, k);
  tuple<intT*, intT*> G;
  t0.start();
  if(directed){
      G = directedKnnGraph(P, A, n, k, sortit);
  }else{
      ;
  }
  cout << "graph-construction-time = " << t0.next() << endl;

  return G;
}

template tuple<intT*, intT*> knnGraph<2>(point<2>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<3>(point<3>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<4>(point<4>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<5>(point<5>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<6>(point<6>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<7>(point<7>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<8>(point<8>*, intT, intT, bool, bool);
template tuple<intT*, intT*> knnGraph<9>(point<9>*, intT, intT, bool, bool);