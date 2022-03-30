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

#define XXH_PRIVATE_API
#include <batchKdtree/shared/geometry.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include <xxHash/xxhash.h>

namespace pargeo::batchKdTree {

#define ATOMIC_BUCKETS
//#define BIT_VECTOR

template <int dim>
class BloomFilter {
  static constexpr int NUM_ARRAYS = 4;
  static constexpr int EMPTY_BUCKETS_GRANULARITY = 1024;
  static constexpr int FILL_BUCKETS_GRANULARITY = 1024;

#ifdef ATOMIC_BUCKETS
#ifdef BIT_VECTOR
  typedef char bucketT;
#else
  typedef bool bucketT;
#endif
  std::atomic<bucketT> *buckets[NUM_ARRAYS];
#else
  std::vector<char> buckets[NUM_ARRAYS];
#endif

  typedef point<dim> pointT;

  const size_t num_buckets;
  const size_t buckets_size;

  size_t hash(const pointT &p, int bucket_num) {
    auto h = (uint64_t)XXH64(p.x, dim * sizeof(double), bucket_num);
    return h % num_buckets;
  }

  static inline bool bit_set(char val, int bit) { return ((val) & (1 << bit)) != 0; }
  static inline char set_bit(char val, int bit) { return ((val) | (1 << bit)); }

  bool get(int j, size_t idx) {
#ifdef ATOMIC_BUCKETS
#ifdef BIT_VECTOR
    return bit_set(buckets[j][idx / 8], idx % 8);
#else
    return buckets[j][idx];
#endif
#else
    return buckets[j][idx];
#endif
  }

  void set(int j, size_t idx) {
    if (!get(j, idx)) {
#ifdef ATOMIC_BUCKETS
#ifdef BIT_VECTOR
      bucketT exp = buckets[j][idx / 8];
      while (!buckets[j][idx / 8].compare_exchange_weak(exp, set_bit(exp, idx % 8)) &&
             !bit_set(exp, idx % 8))
        ;
#else
      bucketT exp = 0;
      buckets[j][idx].compare_exchange_strong(exp, 1);
#endif
#else
      buckets[j][idx] = 1;
#endif
    }
  }

 public:
  BloomFilter(size_t num_points)
      : num_buckets(num_points * 2),
        buckets_size(
#ifdef BIT_VECTOR
            (num_points * 2 + 7) / 8
#else
            num_points * 2
#endif
        ) {
    // TODO: compute better bound for [num_buckets]
    for (int i = 0; i < NUM_ARRAYS; i++) {
#ifdef ATOMIC_BUCKETS
      buckets[i] = new std::atomic<bucketT>[buckets_size];
#else
      buckets[i].resize(buckets_size, 0);
#endif
    }
  }

#ifdef ATOMIC_BUCKETS
  ~BloomFilter() {
    for (int i = 0; i < NUM_ARRAYS; i++) {
      delete[] buckets[i];
    }
  }
#endif

  void clear() {
    // empty buckets
    parlay::parallel_for(0, NUM_ARRAYS, [&](size_t i) {
      parlay::parallel_for(
          0, buckets_size, [&](size_t j) { buckets[i][j] = 0; }, EMPTY_BUCKETS_GRANULARITY);
    });
  }

  template <class R>
  void insert(const R &points) {
    // set buckets
    parlay::parallel_for(
        0,
        points.size(),
        [&](size_t i) {
          for (int j = 0; j < NUM_ARRAYS; j++) {
            auto idx = hash(points[i], j);
            set(j, idx);
          }
        },
        FILL_BUCKETS_GRANULARITY);
  }

  template <class R>
  void build(const R &points) {
    clear();
    insert(points);
  }

  bool might_contain(const point<dim> &p) {
    for (int i = 0; i < NUM_ARRAYS; i++) {
      auto idx = hash(p, i);
      if (!get(i, idx)) return false;
    }
    return true;
  }

  template <class R>
  parlay::sequence<point<dim>> filter(const R &points) {
    return parlay::filter(points, [this](const point<dim> &p) { return this->might_contain(p); });
  }
};

} // End namespace batchKdTree