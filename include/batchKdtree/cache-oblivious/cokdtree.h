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

#include "parlay/parallel.h"
#include "parlay/sequence.h"

#include "utils.h"
#include "../shared/kdnode.h"
#include "../shared/kdtree.h"

#include "../shared/macro.h"

#ifdef ALL_USE_BLOOM
#include "../shared/bloom.h"
#endif

#ifdef PRINT_COKDTREE_TIMINGS
#include "common/get_time.h"
#endif

namespace pargeo::batchKdTree {

// Static Cache-Oblivious Tree
// TODO: add deletes
template <int dim, class objT, bool parallel = false, bool coarsen = false>
class CO_KdTree : public KdTree<dim, objT, parallel, coarsen> {
  typedef KdTree<dim, objT, parallel, coarsen> BaseTree;
  using typename BaseTree::nodeT;

  static const auto leaf_size = (coarsen ? CLUSTER_SIZE : 1);
  // Recursive Build --------------------------------------------------------------------------
  // PARALLEL
  template <bool top>
  void buildKdtRecursiveParallel(parlay::slice<objT *, objT *> items,
                                 nodeT *node_array,
                                 int split_dim,
                                 int num_levels) {
    assert(parallel);
    // DEBUG_MSG("buildKdt" << (top ? "Top" : "Bottom") << "Parallel(items.size() = " <<
    // items.size()
    //<< ", split_dim = " << split_dim << ", num_levels = " << num_levels
    //<< ")");

#ifdef PRINT_COKDTREE_TIMINGS
    bool print_timer = (items.size() > 5000000);
    // bool print_timer = (items.size() > 500000);
#endif

    if (top) {
      // Base case: perform a split
      if (num_levels == 1) {
#ifdef PRINT_COKDTREE_TIMINGS
        timer t("CoKdTree");
        t.start();
#endif
        assert(items.size() > 1);
        auto median = parallelMedianPartition<objT>(items, split_dim);
        assert(node_array[0].isEmpty());
        new (&node_array[0]) nodeT(split_dim, median, items);
#ifdef PRINT_COKDTREE_TIMINGS
        if (print_timer) {
          std::stringstream ss;
          ss << "Placed Split (items.size() = " << items.size() << ")";
          t.report(t.get_next(), ss.str());
        }
#endif
        return;
      }
    } else {
      assert(num_levels > 1);
      // should never hit a base case because of base-case coarsening!
    }

#ifdef PRINT_COKDTREE_TIMINGS
    timer t("CoKdTree");
    t.start();
#endif
    // Recursive case
    int bottom_num_levels = bottomNumLevels(num_levels);
    int top_num_levels = num_levels - bottom_num_levels;
#ifdef PRINT_COKDTREE_TIMINGS
    // print_timer = (top_num_levels == 5) && (bottom_num_levels == 16); // CO-tree 1M
    // print_timer = (top_num_levels == 8) && (bottom_num_levels == 16); // log-tree 10M
    print_timer = (top_num_levels == 9) && (bottom_num_levels == 16);  // CO-tree 10M
#endif

    // Make first call - builds the top tree and partitions [items] on that split.
    auto originalNodeArray = node_array;
    buildKdtTopParallel(items, node_array, split_dim, top_num_levels);

#ifdef PRINT_COKDTREE_TIMINGS
    if (print_timer) {
      std::stringstream ss;
      ss << "Built Top (#top levels, #bottom levels) = (" << top_num_levels << ", "
         << bottom_num_levels << ")";
      t.report(t.get_next(), ss.str());
    }
#endif

    auto top_offset = numNodesTop(top_num_levels);
    node_array += top_offset;

    int num_subtrees = 1 << top_num_levels;
    assert(child_indices[top_num_levels].size() == (size_t)(num_subtrees / 2));
    assert(split_points[top_num_levels].size() == (size_t)num_subtrees);

    size_t size_per_leaf = items.size() / num_subtrees;
    size_t remainder = items.size() % num_subtrees;
    auto left_endpoints = computeLeftEndpoints(top_num_levels, size_per_leaf, remainder);

    parlay::sequence<size_t> nodeOffset;
    if (!top) {
      nodeOffset = computeBottomChildSizes(top_num_levels, size_per_leaf, remainder);
      assert((int)nodeOffset.size() == num_subtrees);
    }
    parlay::parallel_for(
        0,
        num_subtrees,
        [&](size_t i) {
          auto p = i / 2;
          auto parent_idx = child_indices[top_num_levels][p];
          assert(p == i / 2);
          // find endpoints
          auto left_endpoint = left_endpoints[i];
          auto right_endpoint = left_endpoints[i + 1];
          // DEBUG_MSG("size_per_leaf, remainder : " << size_per_leaf << ", " << remainder);
          // DEBUG_MSG("interval: [" << left_endpoint << ", " << right_endpoint << ")");

          auto cur_node_array =
              node_array + (top ? (i * numNodesTop(bottom_num_levels)) : (nodeOffset[i]));
          // assign pointers to this node
          if (i % 2 == 0)
            originalNodeArray[parent_idx].setLeft(cur_node_array);
          else
            originalNodeArray[parent_idx].setRight(cur_node_array);
          // this->parents[cur_node_array - this->nodes] = &originalNodeArray[parent_idx];

          if (top) {
            buildKdtTopParallel(items.cut(left_endpoint, right_endpoint),
                                cur_node_array,
                                (split_dim + top_num_levels) % dim,
                                bottom_num_levels);
          } else {
            buildKdtBottomParallel(items.cut(left_endpoint, right_endpoint),
                                   cur_node_array,
                                   (split_dim + top_num_levels) % dim);
          }
        },
        1);

#ifdef PRINT_COKDTREE_TIMINGS
    if (print_timer) {
      std::stringstream ss;
      ss << "Built Bottom (#top levels, #bottom levels) = (" << top_num_levels << ", "
         << bottom_num_levels << ")";
      t.report(t.get_next(), ss.str());
    }
#endif

#ifndef NDEBUG
    for (int i = 0; i < num_subtrees / 2; i++) {
      auto parent_idx = child_indices[top_num_levels][i];
      auto &parent = originalNodeArray[parent_idx];
      auto left_points = parent.getLeft()->countPoints();
      auto right_points = parent.getRight()->countPoints();
      if (!((left_points == right_points) || (left_points + 1 == right_points))) {
        std::cerr << "ERROR: buildKdt" << (top ? "Top" : "Bottom")
                  << "(items.size() = " << items.size() << ", split_dim = " << split_dim
                  << ", num_levels = " << num_levels << ")Parallel" << std::endl;
        std::cerr << "ERROR: occurs at bottom subtree " << i << " of the above call" << std::endl;
        std::stringstream ss;
        ss << "buildKdt" << (top ? "Top" : "Bottom") << "Parallel!: (left#, right#) = ("
           << left_points << ", " << right_points << ")";
        throw std::runtime_error(ss.str());
      }
    }
#endif
  }

  // Specialized Build Functions
  void buildKdtTopParallel(parlay::slice<objT *, objT *> items,
                           nodeT *node_array,
                           int split_dim,
                           int num_levels) {
    assert(parallel);
    if (buildTopInParallel(num_levels, items.size())) {
      buildKdtRecursiveParallel<true>(items, node_array, split_dim, num_levels);
    } else {
      buildKdtRecursive<true>(items, node_array, split_dim, num_levels);
    }
  }

  void buildKdtBottomParallel(parlay::slice<objT *, objT *> items,
                              nodeT *node_array,
                              int split_dim) {
    assert(parallel);
    int N = items.size();  // number of leaves
    int num_levels = numLevels<coarsen>(N);

    if (buildBottomInParallel(num_levels, items.size())) {
      buildKdtRecursiveParallel<false>(items, node_array, split_dim, num_levels);
    } else {
      buildKdtRecursive<false>(items, node_array, split_dim, num_levels);
    }
  }

  // Recursive Build --------------------------------------------------------------------------
  // SERIAL
  // TODO: figure out if tracking split points/child indices affects cache-obliviousness
  template <bool top>
  size_t buildKdtRecursive(parlay::slice<objT *, objT *> items,
                           nodeT *node_array,
                           int split_dim,
                           int num_levels) {
    // DEBUG_MSG("buildKdt" << (top ? "Top" : "Bottom") << "(items.size() = " << items.size()
    //<< ", split_dim = " << split_dim << ", num_levels = " << num_levels
    //<< ")");

    // Base Cases
    // TODO: consider coarsening the base case (e.g. for bottom:)

    if (top) {
      // Base case: perform a split
      if (num_levels == 1) {
        assert(items.size() > 1);
        auto median = serialMedianPartition<objT>(items, split_dim);
        assert(node_array[0].isEmpty());
        new (&node_array[0]) nodeT(split_dim, median, items);
        return 1;
      }
    } else {
      // Base case: insert the leaves
      if (num_levels == 1) {
        // DEBUG_MSG("CO Leaf: " << items.size() << " points");
        assert(items.size() <= leaf_size);
        assert(items.size() > 0);
        assert(node_array[0].isEmpty());
        // DEBUG_MSG("Placing (" << items[0].coordinate(0) << ", " << items[0].coordinate(1) <<
        // ")");
        new (&node_array[0]) nodeT(items);
        return 1;
      } /*else if (num_levels == 2) {  // coarser base case
        assert(items.size() == 2);
        assert(items[0].coordinate(split_dim) != items[1].coordinate(split_dim));

        auto median = (items[0].coordinate(split_dim) + items[1].coordinate(split_dim)) / 2;
        if (items[0].coordinate(split_dim) > items[1].coordinate(split_dim)) {  // sort
          auto tmp = items[0];
          items[0] = items[1];
          items[1] = tmp;
        }

        // construct the 3-node tree
        auto left = new (&node_array[1]) nodeT(items.cut(0, 1));
        auto right = new (&node_array[2]) nodeT(items.cut(1, 2));
        auto parent = new (&node_array[0]) nodeT(split_dim, median, items);
        parent->setLeft(left);
        parent->setRight(right);
        return 3;
      }*/
    }

    // Recursive case
    int bottom_num_levels = bottomNumLevels(num_levels);
    int top_num_levels = num_levels - bottom_num_levels;

    // Make first call - builds the top tree and partitions [items] on that split.
    auto originalNodeArray = node_array;
    node_array += buildKdtTop(items, node_array, split_dim, top_num_levels);

    int num_subtrees = 1 << top_num_levels;
    assert(child_indices[top_num_levels].size() == (size_t)(num_subtrees / 2));
    assert(split_points[top_num_levels].size() == (size_t)num_subtrees);

#ifndef NDEBUG
    // std::cout << "top splitpoints: [";
    // for (const auto &ep : split_points[top_num_levels])
    // std::cout << ep << ", ";
    // std::cout << "]" << std::endl;
#endif

    size_t left_endpoint = 0;
    size_t size_per_leaf = items.size() / num_subtrees;
    size_t remainder = items.size() % num_subtrees;
    for (int i = 0; i < num_subtrees; i++) {
      // assign pointers to this node
      auto p = i / 2;
      auto parent_idx = child_indices[top_num_levels][p];

      if (i % 2 == 0)
        originalNodeArray[parent_idx].setLeft(node_array);
      else
        originalNodeArray[parent_idx].setRight(node_array);
      // this->parents[node_array - this->nodes] = &originalNodeArray[parent_idx];

      // construct this subtree
      // DEBUG_MSG("size_per_leaf, remainder : " << size_per_leaf << ", " << remainder);
      auto num_in_bucket = size_per_leaf + ((split_points[top_num_levels][i] <= remainder) ? 1 : 0);
      auto right_endpoint = left_endpoint + num_in_bucket;  // exclusive
      // DEBUG_MSG("interval: [" << left_endpoint << ", " << right_endpoint << ")");
      if (top) {
        node_array += buildKdtTop(items.cut(left_endpoint, right_endpoint),
                                  node_array,
                                  (split_dim + top_num_levels) % dim,
                                  bottom_num_levels);
      } else {
        node_array += buildKdtBottom(items.cut(left_endpoint, right_endpoint),
                                     node_array,
                                     (split_dim + top_num_levels) % dim);
      }
      left_endpoint = right_endpoint;
    }

#ifndef NDEBUG
    for (int i = 0; i < num_subtrees / 2; i++) {
      auto parent_idx = child_indices[top_num_levels][i];
      auto &parent = originalNodeArray[parent_idx];
      auto left_points = parent.getLeft()->countPoints();
      auto right_points = parent.getRight()->countPoints();
      if (!((left_points == right_points) || (left_points + 1 == right_points))) {
        std::cerr << "ERROR: buildKdt" << (top ? "Top" : "Bottom")
                  << "(items.size() = " << items.size() << ", split_dim = " << split_dim
                  << ", num_levels = " << num_levels << ")" << std::endl;
        std::cerr << "ERROR: occurs at bottom subtree " << i << " of the above call" << std::endl;
        std::stringstream ss;
        ss << "buildKdt" << (top ? "Top" : "Bottom") << "!: (left#, right#) = (" << left_points
           << ", " << right_points << ")";
        throw std::runtime_error(ss.str());
      }
    }
#endif

    return node_array - originalNodeArray;  // return total number of used spaces
  }

  // Specialized Build Functions
  size_t buildKdtTop(parlay::slice<objT *, objT *> items,
                     nodeT *node_array,
                     int split_dim,
                     int num_levels) {
    return buildKdtRecursive<true>(items, node_array, split_dim, num_levels);
  }

  size_t buildKdtBottom(parlay::slice<objT *, objT *> items, nodeT *node_array, int split_dim) {
    int N = items.size();  // number of leaves
    int num_levels = numLevels<coarsen>(N);

    return buildKdtRecursive<false>(items, node_array, split_dim, num_levels);
  }

  // Base Building Functions
  void buildKdt() {
    // check that it's not degenerate -> at least one split node w/ 2 leaves
    // DEBUG_MSG("buildKdtBottom: items.size() = " << this->items.size()
    //<< ", .size() = " << this->size());
    assert(this->size() > leaf_size);
    assert(numLevels<coarsen>(this->size()) > 1);
    // initialize(numLevels(this->size()));
    if (parallel) {
      buildKdtBottomParallel(this->items.cut(0, this->size()), this->nodes, 0);
    } else {
      buildKdtBottom(this->items.cut(0, this->size()), this->nodes, 0);
    }
  }

  void buildKdt(__attribute__((unused)) parlay::slice<bool *, bool *> flags) {
    throw std::runtime_error(
        "parallel build with flags unimplemented: need some parallel median finding/partition "
        "function!");
  }

 public:
  CO_KdTree(int log2size) : BaseTree(log2size) { initialize(0); }

  // Just a convenience wrapper for tests
  template <class R>
  CO_KdTree(const R &points) : BaseTree(points) {
    initialize(0);
#ifdef PRINT_COKDTREE_TIMINGS
    this->mark_time("Initialize");
#endif
    parlay::sequence<objT> tmp;
    // tmp.resize(points.size());
    // parlay::parallel_for(0, points.size(), [&](size_t i) { tmp[i] = points[i]; });
    tmp.assign(points.begin(), points.end());  // faster than parallel for
#ifdef PRINT_COKDTREE_TIMINGS
    this->mark_time("Point Copy");
#endif
    build(std::move(tmp));
#ifdef PRINT_COKDTREE_TIMINGS
    this->mark_time("Done");
    this->total_time();
#endif
  }

  // MODIFY --------------------------------------------
  /*!
   * Build a new kd-tree in this tree over the input points.
   * This function should only be called on an empty tree.
   * @param points the list of points to build the tree over. Is moved into the tree.
   */
  void build(parlay::sequence<objT> &&points) {
#ifdef PRINT_COKDTREE_TIMINGS
    this->mark_time("Build Called");
#endif
    assert(this->cur_size == 0);
    this->items = points;
    auto build_tree = [&]() {
      auto n = points.size();
      // this assertion holds for build, but not inserts
      // assert(n > this->capacity() / 2);

      // save the input points into [items]
      this->cur_size = n;
      this->build_size = n;
#ifdef PRINT_COKDTREE_TIMINGS
      this->mark_time("Setup");
#endif

      if (n == 0) return;
      buildKdt();
#ifdef PRINT_COKDTREE_TIMINGS
      this->mark_time("Build");
#endif
      this->nodes[0].recomputeBoundingBoxSubtree();  // have to do this afterwards
#ifdef PRINT_COKDTREE_TIMINGS
      this->mark_time("Bounding");
#endif
    };

#ifdef ALL_USE_BLOOM
    auto points_copy = this->items;
    parlay::par_do([&]() { this->bloom_filter.build(points_copy); }, build_tree);
#else
    build_tree();
#endif
    // if (parallel) {
    // auto flags = parlay::sequence<bool>(n);
    // auto flagSlice = parlay::slice(flags.begin(), flags.end());
    // buildKdt(flagSlice);
    //} else {
    // buildKdt();
    //}
  }

  template <class R>
  void insert(const R &points) {
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

    if (this->cur_size == 0) {
      // could be move if points is a sequence
      parlay::sequence<objT> to_insert;
      to_insert.assign(points);
      build(std::move(to_insert));
    } else {
      // gather points from tree
      parlay::sequence<objT> gather(this->cur_size);
      [[maybe_unused]] auto num_moved = this->moveElementsTo(gather.cut(0, this->cur_size));
      assert(num_moved == gather.size());

      // add the new points and rebuild
      gather.append(points);
      build(std::move(gather));
    }
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
    if (rebuild) {
      auto points_copy = points;
      BaseTree::bulk_erase(points_copy);
    } else {
      BaseTree::bulk_erase(points);
    }
#endif

    if (rebuild) {
      parlay::sequence<objT> elements(this->cur_size);
      this->moveElementsTo(elements.cut(0, this->cur_size));
      build(std::move(elements));
    }
  }
};

} // End namespace batchKdTree