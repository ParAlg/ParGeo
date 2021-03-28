#ifndef UNION_FIND_H
#define UNION_FIND_H

//adapted from PBBS
// #include "common/pbbs/utils.h"
#include <limits>
#include <iostream>
using namespace std;

typedef int intT;
typedef unsigned int uintT;
typedef double floatT;
static intT intMax() {return numeric_limits<intT>::max();}
static uintT uintMax() {return numeric_limits<uintT>::max();}
static floatT floatMax() {return numeric_limits<floatT>::max();}
static floatT floatMin() {return numeric_limits<floatT>::lowest();}
#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))
// #include "parallel.h"
#include "parlay/parallel.h" // can't use the unionFind in common, parallel conflict

struct unionFind {

  int *parents;
  int *hooks;

  unionFind(int n) {
    parents = newA(int, n);
    parlay::parallel_for (0, n, [&](int i) {parents[i] = intMax();});
    hooks = newA(int, n);
    parlay::parallel_for (0, n, [&](int i) {hooks[i] = intMax();});
  }

  void del() {free(parents);}

  inline int find(int i) {
    int j = i;
    if (parents[j] == intMax()) return j;
    do j = parents[j];
    while (parents[j] < intMax());
    int tmp;
    while((tmp=parents[i])<j){ parents[i]=j; i=tmp;}
    return j;
  }

  void link(int u, int v) {
    while(1){
      u = find(u);
      v = find(v);
      if(u == v) break;
      if(u > v) swap(u,v);
      if(hooks[u] == intMax() && __sync_bool_compare_and_swap(&hooks[u], intMax(), u)){
        parents[u]=v;
        break;
      }}
  }
};

struct edgeUnionFind {
  typedef pair<int, int> edgeT;
  int *parents;
  edgeT *hooks;
  int n;

edgeUnionFind(int nn): n(nn) {
  parents = newA(int, n);
  parlay::parallel_for (0, n, [&](int i) {parents[i] = intMax();});
  hooks = newA(edgeT, n);
  parlay::parallel_for (0, n, [&](int i) {
      hooks[i] = make_pair(intMax(), intMax());});
}

  void del() {free(parents);}

  inline int find(int i) {
    int j = i;
    if (parents[j] == intMax()) return j;
    do j = parents[j];
    while (parents[j] < intMax());
    int tmp;
    while((tmp=parents[i])<j) {parents[i]=j; i=tmp;} 
    return j;
  }

  int link(int u, int v) {
    int c_from = u;
    int c_to = v;
    while(1) {
      u = find(u);
      v = find(v);
      if(u == v) break;
      if(u > v) swap(u,v);
      if(hooks[u].first == intMax() && __sync_bool_compare_and_swap(&hooks[u].first, intMax(), c_from)){
        parents[u]=v;
        hooks[u].second=c_to;
        break;
      }
    }
    return parents[u];
  }

  edgeT getEdge(int idx) {
    return hooks[idx];
  }

};
  

#endif

