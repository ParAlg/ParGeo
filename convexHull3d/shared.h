#pragma once

#include "parlay/hash_table.h"

// Example for hashing numeric values.
// T must be some integer type
template <class T>
struct hash_pointer {
  using eType = T;
  using kType = T;
  eType empty() { return nullptr; }
  kType getKey(eType v) { return v; }
  size_t hash(kType v) { return static_cast<size_t>(hash64(size_t(v))); }
  int cmp(kType v, kType b) { return (v > b) ? 1 : ((v == b) ? 0 : -1); }
  bool replaceQ(eType, eType) { return 0; }
  eType update(eType v, eType) { return v; }
  bool cas(eType* p, eType o, eType n) {
    // TODO: Make this use atomics properly. This is a quick
    // fix to get around the fact that the hashtable does
    // not use atomics. This will break for types that
    // do not inline perfectly inside a std::atomic (i.e.,
    // any type that the standard library chooses to lock)
    return std::atomic_compare_exchange_strong_explicit(
      reinterpret_cast<std::atomic<eType>*>(p), &o, n, std::memory_order_relaxed, std::memory_order_relaxed);
  }
};
