#ifndef KDTREE_SHARED_KNNBUFFER_H
#define KDTREE_SHARED_KNNBUFFER_H

// TAKEN with slight modifications FROM
// https://github.mit.edu/yiqiuw/pargeo/blob/master/knnSearch/kdTree/kdtKnn.h Later, need to merge +
// refer to that rather than copying here
// #include <common/geometry.h>
#include "pargeo/point.h"
namespace knnBuf {

typedef int intT;
typedef double floatT;

template <typename T>
struct elem {
  floatT cost;  // Non-negative
  T entry;
  elem(floatT t_cost, T t_entry) : cost(t_cost), entry(t_entry) {}
  elem() : cost(std::numeric_limits<floatT>::max()) {}
  bool operator<(const elem& b) const {
    if (cost < b.cost) return true;
    return false;
  }
};

template <typename T>
struct buffer {
  typedef parlay::slice<elem<T>*, elem<T>*> sliceT;
  /*const*/ intT k;  // not const because of assignment in dualKnn
  intT ptr;
  sliceT buf;
  bool cached;

  buffer() : k(-1), ptr(0), buf(), cached(false) {}
  buffer(intT t_k, sliceT t_buf) : k(t_k), ptr(0), buf(t_buf), cached(true) {}

  inline void reset() { ptr = 0; }

  bool hasK() { return ptr >= k; }

  elem<T> keepK() {
    if (ptr < k) throw std::runtime_error("Error, kbuffer not enough k.");
    if (!cached) {  // only need to do this if modified since last call
      std::nth_element(buf.begin(), buf.begin() + k - 1, buf.begin() + ptr);
      ptr = k;
      cached = true;
    }
    return buf[k - 1];
  }

  void insert(elem<T> t_elem) {
    buf[ptr++] = t_elem;
    cached = false;
    if (ptr >= (int)buf.size()) keepK();
  }

  elem<T> operator[](intT i) {
    if (i < ptr)
      return buf[i];
    else
      return elem<T>();
  }
};

template <int dim>
parlay::sequence<const pargeo::point<dim>*> bruteforceKnn(const parlay::sequence<pargeo::point<dim>>& queries,
                                                  size_t k) {
  auto out = parlay::sequence<elem<const pargeo::point<dim>*>>(2 * k * queries.size());
  auto idx = parlay::sequence<const pargeo::point<dim>*>(k * queries.size());
  parlay::parallel_for(0, queries.size(), [&](size_t i) {
    auto q = queries[i];
    buffer buf = buffer<const pargeo::point<dim>*>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
    for (intT j = 0; j < (int)queries.size(); ++j) {
      auto p = &queries[j];
      buf.insert(elem(q.dist(p), p));
    }
    buf.keepK();

    // store results
    for (size_t j = 0; j < k; j++) {
      idx[i * k + j] = buf[j].entry;
    }
  });
  return idx;
}

}  // namespace knnBuf

#endif  //  KDTREE_SHARED_KNNBUFFER_H
