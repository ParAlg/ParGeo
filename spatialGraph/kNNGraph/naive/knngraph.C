#include "kNearestNeighbor/kdTree/knn.C"
#include "pbbs/gettime.h"
#include "pbbs/sampleSort.h"
#include "pbbs/utils.h"
#include "pbbs/sequence.h"


#include <cstdint>
using namespace std;

// #define DEBUG
#define getNN(i, kold, j) (reinterpret_cast<uintptr_t>(A[i*kold+j]) - base_addr) / point_size // the jth nn of point i
#define distNN(i, j) P[i].pointDist(*A[i*kold+j]) //dist from point i to its jth nn

void printDirectedGraph(intT* offsets, intT* edges, intT n, intT k){
#ifdef DEBUG
    // for(intT i = 0; i < n; ++i){
    //     cout << offsets[i] << endl;
    // }
    // cout << "======" << endl;
    for(intT i = 0; i < n; ++i){
        cout  << ":::" <<  i << ":::" << endl;
        for(intT j = 0; j < k; ++j){
            cout << edges[i*k + j] << endl;
        }
    }
    cout << "======" << endl;
#endif
}

void printUndirectedGraph(intT* offsets, intT* edges, intT n, intT m){
#ifdef DEBUG
    for(intT i = 0; i < n; ++i){
        cout << offsets[i] << endl;
    }
    cout << "======" << endl;
    for(intT i = 0; i < m; ++i){
        cout << edges[i] << endl;
    }
    cout << "======" << endl;
#endif
}

// A: chunks of k NNs
// doo not include self edge
template<int dim>
tuple<intT*, intT*> directedKnnGraph(point<dim>* P, point<dim>** A, intT n, intT k, bool sortit){
    typedef point<dim> pointT;
    intT *offsets = (intT*) malloc((n)*sizeof(intT));
    intT *edges = (intT*) malloc((n*k)*sizeof(intT));
    auto base_addr = reinterpret_cast<uintptr_t>(P);
    auto point_size = sizeof(pointT);
    k = k-1;
    parallel_for (0, n, [&](intT i) {offsets[i] = k;});
    sequence::prefixSum(offsets, 0, n);
    parallel_for (0, n, [&](intT i) {
        parallel_for (1, k+1, [&](intT j) {
            edges[i*k + j-1] = (reinterpret_cast<uintptr_t>(A[i*(k+1)+j]) - base_addr) / point_size;
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

    printDirectedGraph(offsets, edges, n, k);

    return make_tuple(offsets, edges);
}



// #define CACHELINE_SIZE 128
//upper bound on other points in its knn?
// if equal distance, does knn graph include the edges?
// if distance equal, search back
template<int dim>
tuple<intT*, intT*> undirectedKnnGraph(point<dim>* P, point<dim>** A, intT n, intT k, bool sortit){
    typedef point<dim> pointT;
    auto base_addr = reinterpret_cast<uintptr_t>(P);
    auto point_size = sizeof(pointT);

    intT kold = k;
    k = k-1;
    intT *offsets = (intT*) malloc((n)*sizeof(intT));
    parallel_for (0, n, [&](intT i) {offsets[i] = k;});
    // at location ind = i*k+j-1, -1 if do not need to write the jth nn of point i twice
    // if not -1, write the jth nn of point i to be nn's edges_offsets[ind] neighbor
    intT *edges_offsets = (intT*) malloc((2*n*k)*sizeof(intT)); 


    double *radius = (double *) malloc((n)*sizeof(double));
    parallel_for (0, n, [&](intT i) {
        radius[i] = distNN(i, k); //P[i].pointDist(*A[(i+1)*kold-1]);
    });

    parallel_for (0, n, [&](intT i) {
        parallel_for (1, kold, [&](intT j) {
//    for(intT i = 0; i < n; ++i){
//         for(intT j = 1; j < kold; ++j){
            intT nn = getNN(i, kold, j);
            intT edges_offsets_ind = -1; //default: i in nn's knn
            double dist = distNN(i, j);
            if(dist > radius[nn]){ // i not in nn's knn
                edges_offsets_ind = utils::fetchAndAdd(&offsets[nn], 1);
            }else if(dist == radius[nn]){ //i may or may not be in nn's knn
                intT check_idx = k;
                intT nnn = getNN(nn, kold, check_idx);
                //can optimize if assume A is sorted by id
                while(check_idx > 0 && nnn != i && dist == distNN(nn, check_idx)){ 
                    check_idx --;
                    nnn = getNN(nn, kold, check_idx);
                }
                if(nnn != i){ // i not in nn's knn
                    edges_offsets_ind = utils::fetchAndAdd(&offsets[nn], 1);
                }
            }
            edges_offsets[i*k+j-1] = edges_offsets_ind;
        });
    });
    //     }
    // }

    intT m = sequence::prefixSum(offsets, 0, n);
    intT *edges = (intT*) malloc(m*sizeof(intT));

    parallel_for (0, n, [&](intT i) {
        parallel_for (1, kold, [&](intT j) {
//    for(intT i = 0; i < n; ++i){
//         for(intT j = 1; j < kold; ++j){
            intT nn = getNN(i, kold, j);
            edges[offsets[i] + j-1] = nn;
            intT ind = edges_offsets[i*k+j-1];
            if(ind != -1){ // i not in nn's knn
                edges[offsets[nn] + ind] = i;
            }
        });
    });
    //     }
    // }

    free(edges_offsets);
    free(radius);
    
    printUndirectedGraph(offsets, edges, n, m);

    return make_tuple(offsets, edges);
}
// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
tuple<intT*, intT*> knnGraph(point<dim>* P, intT n, intT k, bool directed, bool sortit) {
  typedef point<dim> pointT;
  if (n < 2) abort();

  timing t0;
  pointT** A = knn(P, n, k);
  tuple<intT*, intT*> G;
  t0.start();
  if(directed){
      cout << k << "-nearest neighbor directed graph construction, " << n << ", dim " << dim << " points" << endl;
      G = directedKnnGraph(P, A, n, k, sortit);
  }else{
      cout << k << "-nearest neighbor undirected graph construction, " << n << ", dim " << dim << " points" << endl;
      G = undirectedKnnGraph(P, A, n, k, sortit);
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