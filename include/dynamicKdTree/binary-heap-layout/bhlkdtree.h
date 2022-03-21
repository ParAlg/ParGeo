#ifndef BHLKDTREE_H
#define BHLKDTREE_H

#include "parlay/parallel.h"
#include "parlay/sequence.h"

#include "../shared/kdnode.h"
#include "../shared/kdtree.h"
#include "../shared/utils.h"

#include "../shared/macro.h"

inline bool buildInParallel(size_t num_points) { return num_points >= BHL_BUILD_BASE_CASE; }

// Dynamic (insert + delete) kd-tree with binary-heap layout
template <int dim, class objT, bool parallel = false, bool coarsen = false>
class BHL_KdTree : public KdTree<dim, objT, parallel, coarsen> {
  typedef KdTree<dim, objT, parallel, coarsen> BaseTree;
  using typename BaseTree::nodeT;

  static const auto leaf_size = (coarsen ? CLUSTER_SIZE : 1);

#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
  void buildKdtRecursive(parlay::slice<objT *, objT *> items, int node_idx, int split_dim)
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
  void buildKdtRecursive(parlay::slice<objT *, objT *> items,
                         parlay::slice<bool *, bool *> flags,
                         int node_idx,
                         int split_dim)
#endif
  {
    if (items.size() <= leaf_size) {  // Base Case
      assert(items.size() > 0);
      new (&this->nodes[node_idx]) nodeT(items);
      return;
    }

    // Recursive Case
    bool parallelBuild = parallel && buildInParallel(items.size());  // should we parallelize?

    // Make the split
    double median;
    size_t right_start;

#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
    if (parallelBuild) {
      median = parallelMedianPartition<objT>(items, split_dim);
    } else {
      median = serialMedianPartition<objT>(items, split_dim);
    }
    right_start = items.size() / 2;
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
    if (parallelBuild) {
      median = parallelSpatialPartition<objT>(items, flags, split_dim, right_start);
    } else {
      median = serialSpatialPartition<objT>(items, split_dim, right_start);
    }
#else
    throw std::runtime_error("invalid partition type");
#endif

    assert(right_start > 0);
    // if (right_start == items.size()) {
    // std::cout << "(split_dim, median, items.size()) = (" << split_dim << ", " << median << ", "
    //<< items.size() << ")" << std::endl;
    //}
    assert(right_start < items.size());

    // place the split
    assert(node_idx >= 0);
    assert((size_t)node_idx < this->num_nodes());
    assert(this->nodes[node_idx].isEmpty());
    auto p = new (&this->nodes[node_idx]) nodeT(split_dim, median, items);

    // Child functors
    auto next_split_dim = (split_dim + 1) % dim;
    auto left_f = [&]() {
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
      auto left_idx = 2 * node_idx + 1;
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
      auto left_idx = node_idx + 1;
#endif
      p->setLeft(&this->nodes[left_idx]);
      // this->parents[left_idx] = p;
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
      buildKdtRecursive(items.cut(0, right_start), left_idx, next_split_dim);
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
      buildKdtRecursive(
          items.cut(0, right_start), flags.cut(0, right_start), left_idx, next_split_dim);
#endif
    };
    auto right_f = [&]() {
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
      auto right_idx = 2 * node_idx + 2;
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
      auto right_idx = node_idx + (2 * right_start);
#endif
      p->setRight(&this->nodes[right_idx]);
      // this->parents[right_idx] = p;
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
      buildKdtRecursive(items.cut(right_start, items.size()), right_idx, next_split_dim);
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
      buildKdtRecursive(items.cut(right_start, items.size()),
                        flags.cut(right_start, flags.size()),
                        right_idx,
                        next_split_dim);
#endif
    };

    // Build Children
    if (parallelBuild) {
      parlay::par_do(left_f, right_f);
    } else {
      left_f();
      right_f();
    }
    assert(!p->isLeaf());
    p->recomputeBoundingBox();
  }

  // Base Building Functions
  void buildKdt() {
    assert(this->size() > leaf_size);  // so it's not degenerate
    buildKdtRecursive(parlay::slice(this->items.begin(), this->items.begin() + this->size()), 0, 0);
  }

  void buildKdt(parlay::slice<bool *, bool *> flags) {
    assert(this->size() > leaf_size);  // so it's not degenerate
    buildKdtRecursive(
        parlay::slice(this->items.begin(), this->items.begin() + this->size()), flags, 0, 0);
  }

  // cur_size, build_size, max_size, items are set before this is called
  void build() {
    assert(this->cur_size != 0);
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
    buildKdt();
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
    auto n = this->items.size();
    auto flags = parlay::sequence<bool>(n);
    auto flagSlice = parlay::slice(flags.begin(), flags.end());
    buildKdt(flagSlice);
#endif
  }

 public:
  BHL_KdTree(int log2size, bool initialize = true) : BaseTree(log2size) {
    if (initialize) this->items = parlay::sequence<objT>(this->max_size);
  }
  BHL_KdTree(const parlay::slice<const objT *, const objT *> &points) : BaseTree(points) {
    this->items = parlay::sequence<objT>(this->max_size);
    build(points);
  }
  BHL_KdTree(const parlay::sequence<objT> &points) : BHL_KdTree(points.cut(0, points.size())) {}

  // MODIFY --------------------------------------------
  /*!
   * Build a new kd-tree in this tree over the input points.
   * This function should only be called on an empty tree.
   * @param points the list of points to build the tree over.
   */
  void build_no_bloom(const parlay::slice<const objT *, const objT *> &points) {
    assert(this->cur_size == 0);
    size_t n = points.size();

    // save the input points into [items]
    // TODO: does this parallelize? items.assign(points.begin(), points.end());
    this->cur_size = n;
    this->build_size = n;
    parlay::parallel_for(0, n, [&](size_t i) { this->items[i] = points[i]; });

    build();
  }

  void build(const parlay::slice<const objT *, const objT *> &points) {
    assert(this->cur_size == 0);

#ifdef ALL_USE_BLOOM
    parlay::par_do([&]() { this->bloom_filter.build(points); }, [&]() { build_no_bloom(points); });
#else
    build_no_bloom(points);
#endif
  }

  // should somehow indicate that inserts aren't allowed
  void build(parlay::sequence<objT> &&points) {
    assert(this->cur_size == 0);
    size_t n = points.size();
    this->cur_size = n;
    this->build_size = n;
    this->items = points;
#ifdef ALL_USE_BLOOM
    auto points_copy = this->items;
    parlay::par_do([&]() { this->bloom_filter.build(points_copy); }, [&]() { build(); });
#else
    build();
#endif
  }

  void insert(const parlay::slice<const objT *, const objT *> &points) {
    if (points.size() == 0) return;
#ifndef NDEBUG
#warning "this is probably bad for runtime - remove later"
    if (points.size() > this->capacity() - this->cur_size) {
      std::stringstream ss;
      ss << "Invalid Insert: (points.size(), max_size, cur_size) = (" << points.size() << ", "
         << this->capacity() << ", " << this->cur_size << ")";
      throw std::runtime_error(ss.str());
    }
#endif

    auto insert_no_bloom = [&]() {
      if (this->cur_size == 0) {
        build_no_bloom(points);
      } else {
        // gather points from tree
        parlay::sequence<objT> gather(this->cur_size);
        [[maybe_unused]] auto num_moved = this->moveElementsTo(gather.cut(0, this->cur_size));
        assert(num_moved == gather.size());

        parlay::parallel_for(0, gather.size() + points.size(), [&](size_t i) {
          if (i < gather.size())
            this->items[i] = gather[i];  // reinsert old items
          else
            this->items[i] = points[i - gather.size()];  // insert new items
        });

        // set sizes
        this->cur_size = gather.size() + points.size();
        this->build_size = this->cur_size;
        build();
      }
    };

#ifdef ALL_USE_BLOOM
    parlay::par_do([&]() { this->bloom_filter.insert(points); }, [&]() { insert_no_bloom(); });
#else
    insert_no_bloom();
#endif
  }

  template <bool rebuild = true>
#ifdef ALL_USE_BLOOM
  void bulk_erase(const parlay::sequence<objT> &points)
#else
  void bulk_erase(parlay::sequence<objT> &points)
#endif
  {
#ifdef ALL_USE_BLOOM
    BaseTree::bulk_erase(points);
#else
    if (rebuild) {  // NOT logtree!
      auto points_copy = points;
      BaseTree::bulk_erase(points_copy);
    } else {
      BaseTree::bulk_erase(points);
    }
#endif

    if (rebuild) {
      auto cursize = this->cur_size;
      parlay::sequence<objT> elements(cursize);
      this->moveElementsTo(elements.cut(0, cursize));
      const auto &const_elements = elements;
      build(const_elements.cut(0, cursize));
    }
  }
};

#endif  // BHLKDTREE_H
