// This code is part of the project "Parallel Batch-Dynamic Kd-Trees"
// Copyright (c) 2021-2022 Rahul Yesantharao, Yiqiu Wang, Laxman Dhulipala, Julian Shun
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

// TAKEN with slight modifications FROM
// https://github.mit.edu/yiqiuw/pargeo/blob/master/knnSearch/kdTree/kdtKnn.h Later, need to merge +
// refer to that rather than copying here
#include "batchKdtree/shared/geometry.h"
namespace pargeo::batchKdTree::knnBuf {

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
parlay::sequence<const point<dim>*> bruteforceKnn(const parlay::sequence<point<dim>>& queries,
                                                  size_t k) {
  auto out = parlay::sequence<elem<const point<dim>*>>(2 * k * queries.size());
  auto idx = parlay::sequence<const point<dim>*>(k * queries.size());
  parlay::parallel_for(0, queries.size(), [&](size_t i) {
    auto q = queries[i];
    buffer buf = buffer<const point<dim>*>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
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

}  // namespace batchKdTree::knnBuf
