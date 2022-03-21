#ifndef KDTREE_SHARED_UTILS_H
#define KDTREE_SHARED_UTILS_H

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"
#include "parlay/monoid.h"
#include "parlay/delayed_sequence.h"

#include "macro.h"

#ifdef DEBUG
#define DEBUG_MSG(str)             \
  do {                             \
    std::cout << str << std::endl; \
  } while (false)
#else
#define DEBUG_MSG(str) \
  do {                 \
  } while (false)
#endif

// MEDIAN FINDING FUNCTIONS ------------------------------------------------------------------------
// convention: always puts more on the right hand side (if odd number of items)
// TODO: think about cache-obliviousness?

template <class objT>
static inline int split_n(const parlay::slice<objT *, objT *> &items) {
  return items.size() / 2;
}
template <class objT, bool select = false>
static inline double split_val(const parlay::slice<objT *, objT *> &items, int dimension) {
  auto split_point = split_n(items);
  if (select) {  // only have split_point in place, not element before, so can only use it.
    return items[split_point].coordinate(dimension);
  } else {
    return (items[split_point - 1].coordinate(dimension) +
            items[split_point].coordinate(dimension)) /
           2.0;
  }
}

// Implementations: rearrange [items] so that the median is in the middle ---------------------
// sort based implementation
template <class objT, bool parallel>
static inline void medianPartitionSort(parlay::slice<objT *, objT *> &items, int dimension) {
  auto compare = [dimension](const objT &l, const objT &r) {
    return l.coordinate(dimension) < r.coordinate(dimension);
  };

  if (parallel) {
    parlay::sort_inplace(items, compare);
  } else {
    std::sort(items.begin(), items.end(), compare);
  }
}

// serial select based implementation
template <class objT>
static inline void serialMedianPartitionSelect(parlay::slice<objT *, objT *> &items,
                                               int dimension) {
  auto compare = [dimension](const objT &l, const objT &r) {
    return l.coordinate(dimension) < r.coordinate(dimension);
  };

  auto split_point = split_n(items);
  std::nth_element(items.begin(), items.begin() + split_point, items.end(), compare);
}

// Top level functions ------------------
// object median
template <class objT>
double parallelMedianPartition(parlay::slice<objT *, objT *> items, int dimension) {
#ifdef USE_MEDIAN_SELECTION
  serialMedianPartitionSelect<objT>(items, dimension);  // selection-based version
  return split_val<objT, true>(items, dimension);
#else
  medianPartitionSort<objT, true>(items, dimension);  // sort-based version
  return split_val(items, dimension);
#endif
}

template <class objT>
double serialMedianPartition(parlay::slice<objT *, objT *> items, int dimension) {
  //#ifdef USE_MEDIAN_SELECTION
  serialMedianPartitionSelect<objT>(items, dimension);  // selection-based version
  return split_val<objT, true>(items, dimension);
  //#else
  // medianPartitionSort<objT, false>(items, dimension);  // sort-based version
  // return split_val(items, dimension);
  //#endif
}

// PARTITION FUNCTIONS -----------------------------------------------------------------------------

template <class objT>
auto serialPartition(parlay::slice<objT *, objT *> points, int dimension, double value) {
  auto comp_lambda = [dimension, value](const objT &o) { return o.coordinate(dimension) < value; };
  auto ret_it = std::partition(points.begin(), points.end(), comp_lambda);
  return ret_it - points.begin();
}

template <class objT>
auto parallelPartition(parlay::slice<objT *, objT *> points,
                       parlay::slice<bool *, bool *> flags,
                       int dimension,
                       double value) {
  // construct partition flags (split_two puts false first)
#ifdef PRINT_PARALLEL_PARTITION_TIMINGS
  bool print_timings = (points.size() == 2000000);
  timer t("parallelPartition");
  auto mtime = [&](const std::string &s) {
    if (print_timings) t.report(t.get_next(), s);
  };
#endif

  parlay::parallel_for(0, points.size(), [dimension, value, points, flags](size_t i) {
    if (points[i].coordinate(dimension) < value)
      flags[i] = false;
    else
      flags[i] = true;
  });

#ifdef PRINT_PARALLEL_PARTITION_TIMINGS
  mtime("flags set");
#endif

  // perform partition and copy back to original points array
  const auto &[res, split_pt] = parlay::internal::split_two(points, flags);

#ifdef PRINT_PARALLEL_PARTITION_TIMINGS
  mtime("split2: " + std::to_string(res.size()));
#endif

  const size_t num_blocks = parlay::num_workers() * 8;
  const size_t block_size = (points.size() + num_blocks - 1) / num_blocks;
  parlay::parallel_for(0, num_blocks, [&](size_t i) {
    auto start_idx = i * block_size;
    auto end_idx = std::min((i + 1) * block_size, points.size());
    for (size_t j = start_idx; j < end_idx; j++)
      points[j] = res[j];
  });

#ifdef PRINT_PARALLEL_PARTITION_TIMINGS
  mtime("copyback");
#endif

  return split_pt;
}

bool eraseInParallel(size_t num_points) { return num_points >= ERASE_BASE_CASE; }

template <class TT>
struct minmaxm {
  using T = std::pair<TT, TT>;
  minmaxm() : identity(T(parlay::highest<TT>(), parlay::lowest<TT>())) {}
  T identity;
  static T f(T a, T b) { return T((std::min)(a.first, b.first), (std::max)(a.second, b.second)); }
};

// other partition functions
template <class objT>
double parallelSpatialPartition(parlay::slice<objT *, objT *> items,
                                parlay::slice<bool *, bool *> flags,
                                int dimension,
                                size_t &split_pt) {
  // compute median
  // construct sequence of just the relevant coordinates
  auto f = parlay::delayed_seq<std::pair<double, double>>(items.size(), [&](size_t i) {
    return std::make_pair<double, double>(items[i].coordinate(dimension),
                                          items[i].coordinate(dimension));
  });
  // find min, max - could do this without delayed_seq if we had a full reduce func (different
  // element, return type)
  auto [imin, imax] = parlay::reduce(f, minmaxm<double>());
  auto imed = (imin + imax) / 2;

  // partition
  split_pt = parallelPartition(items, flags, dimension, imed);
  return imed;
}

template <class objT>
double serialSpatialPartition(parlay::slice<objT *, objT *> items,
                              int dimension,
                              size_t &split_pt) {
  // compute median
  double imin = std::numeric_limits<double>::max();
  double imax = std::numeric_limits<double>::lowest();
  assert(imax < -2000);
  for (const auto &pt : items) {
    double val = pt.coordinate(dimension);
    imin = std::min(imin, val);
    imax = std::max(imax, val);
  }
  auto imed = (imin + imax) / 2;

  // partition
  split_pt = serialPartition(items, dimension, imed);
  // if (split_pt == items.size()) {
  // std::cout << "PROBLEM" << std::endl;
  // std::cout << "(split_dim, imin, imax, imed, items.size()) = (" << dimension << ", " << imin
  //<< ", " << imax << ", " << imed << ", " << items.size() << ")" << std::endl;
  // std::vector<double> vals;
  // for (const auto &pt : items) {
  // vals.push_back(pt.coordinate(dimension));
  //}
  // std::sort(vals.begin(), vals.end());
  // for (const auto &val : vals) {
  // std::cout << val << ", ";
  //}
  // std::cout << std::endl;
  // assert(false);
  //}
  return imed;
}

#endif  // KDTREE_SHARED_UTILS_H
