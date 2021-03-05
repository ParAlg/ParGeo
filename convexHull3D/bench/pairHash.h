#pragma once

#include "common/sparse_table.h"

using namespace std;

struct emptyT {};

struct hashIntPair {
  inline size_t operator () (const pair<int, int>& t) const {
    size_t l = t.first;
    size_t r = t.second;
    size_t key = (l << 32) + r;
    return parlay::hash64_2(key);
  }
};

struct pairHash {
  typedef pbbs::sparse_table<pair<int,int>, emptyT, hashIntPair> tableT;

  tableT edgeTable;
  pairHash(size_t n)
    : edgeTable(tableT(n, make_tuple(make_pair<int,int>(-1,-1), emptyT()), hashIntPair())) {
  }

  void mark(int p1, int p2) {
    edgeTable.insert(make_tuple(pair(p1,p2),emptyT()));
  };

  bool processed(int p1, int p2) {
    return edgeTable.contains(pair(p1,p2));
  };
};
