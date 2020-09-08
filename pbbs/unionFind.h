//adapted from pbbs code
#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <limits>

struct edgeUnionFind {
  typedef pair<intT, intT> edgeT;
  intT *parents;
  //intT *hooks; // stores u
  //intT *connect_to; // stores v
  edgeT *hooks;
  intT n;

  // initialize with all roots marked with -1
  edgeUnionFind(intT nn): n(nn) {
    parents = newA(intT, n);
    par_for (intT i=0; i < n; i++) parents[i] = intMax();
    //hooks = (edgeT *) malloc(sizeof(edgeT) * n);
    hooks = newA(edgeT, n);
    par_for (intT i=0; i < n; i++) {
      hooks[i] = make_pair(intMax(), intMax());}
    // hooks = newA(intT, n);
    // parallel_for (intT i=0; i < n; i++) hooks[i] = std::numeric_limits<intT>::max();
    // connect_to = newA(intT, n);
    // parallel_for (intT i=0; i < n; i++) connect_to[i] = std::numeric_limits<intT>::max();
  }

  void del() {free(parents);}

  // Assumes root is negative 
  // Not making parent array volatile improves
  // performance and doesn't affect correctness
  inline intT find(intT i) {
    intT j = i;
    if (parents[j] == intMax()) return j;
    do j = parents[j];
    while (parents[j] < intMax());
    //note: path compression can happen in parallel in the same tree, so
    //only link from smaller to larger to avoid cycles
    intT tmp;
    while((tmp=parents[i])<j) {parents[i]=j; i=tmp;} 
    return j;
  }

  intT link(intT u, intT v) {
    intT c_from = u;
    intT c_to = v;
    while(1) {
      u = find(u);
      v = find(v);
      if(u == v) break;
      if(u > v) swap(u,v);
      //if successful, store the ID of the edge used in hooks[u]
      if(hooks[u].first == intMax() && __sync_bool_compare_and_swap(&hooks[u].first, intMax(), c_from)){
        parents[u]=v;
        hooks[u].second=c_to;
        //connect_to[u]=c_to;
        break;
      }
    }
    return parents[u];
  }

  edgeT getEdge(intT idx) {
    //return make_pair(hooks[idx], connect_to[idx]);
    return hooks[idx];
  }

};
  

#endif

