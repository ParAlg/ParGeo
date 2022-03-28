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
#include "parlay/primitives.h"
#include "../shared/utils.h"

namespace pargeo::batchKdTree {

void funnelSort() {
  // TODO: cache-oblivious sorting
}

// kd-tree structure
template <bool coarsen>
inline int numLevels(int num_leaves) {
  if (coarsen) num_leaves = (num_leaves + (CLUSTER_SIZE - 1)) / CLUSTER_SIZE;
  return 1 + (int)std::ceil(std::log2(num_leaves));  // number of levels
}
inline int numNodesTop(int num_levels) { return (1 << num_levels) - 1; }
inline int numNodesBottom(size_t num_points) { return 2 * num_points - 1; }

inline bool buildTopInParallel(__attribute__((unused)) int num_levels, size_t num_points) {
  if (num_points < CO_TOP_BUILD_BASE_CASE) return false;
  return true;
}
inline bool buildBottomInParallel(__attribute__((unused)) int num_levels, size_t num_points) {
  if (num_points < CO_BOTTOM_BUILD_BASE_CASE) return false;
  return true;
}

inline int hyperCeiling(int n) { return 1 << (int)std::ceil(std::log2(n)); }
inline int bottomNumLevels(int total_num_levels) {
  return hyperCeiling((total_num_levels + 1) / 2);
}

// compute split points in closed form; see memoized version below
std::vector<size_t> computeSplitPoints(int num_levels, int num_leaves) {
  // roughly: the split_points[level] is the order in which elements are placed into this level
  auto num_bottom = (1 << (num_levels));
  std::vector<size_t> split_points(num_bottom + 1);
  split_points[0] = 0;
  for (int level = 0; level < num_levels; level++) {
    auto offset = (1 << (num_levels - 1 - level)) - 1;
    auto adder = (1 << (num_levels - level));
    for (int node = 0; node < (1 << level); node++) {
      auto cur_index = adder * node + offset + 1;
      auto left_endpoint = (node == 0) ? 0 : split_points[cur_index - (offset + 1)];
      auto right_endpoint =
          (node == (1 << level) - 1) ? num_leaves : split_points[cur_index + (offset + 1)];
      split_points[cur_index] = (left_endpoint + right_endpoint) / 2;
    }
  }
  split_points[num_bottom] = num_leaves;
  return split_points;
}

std::vector<std::vector<size_t>> split_points;
void initializeSplitPoints(int num_levels) {
  assert(num_levels > 1);
  if (num_levels < (int)split_points.size()) return;  // already done

  int level = (int)split_points.size();  // the next level to start populating (because of the +1)
  split_points.resize(num_levels + 1);   // allocate space for the levels

  // base case
  if (level == 0) {
    split_points[0] = {1};
    level++;
  }

  // Note: this can be computed in closed-form, but it seems like memoizing makes sense?
  for (; level <= num_levels; level++) {
    auto num_leaves = (1 << level);
    split_points[level].resize(num_leaves);
    for (int leaf = 0; leaf < num_leaves; leaf++) {
      split_points[level][leaf] =
          split_points[level - 1][leaf / 2] + ((leaf % 2 == 0) ? (1 << (level - 1)) : 0);
    }
  }
}

// compute left indices of children at a given level with a given remainder
parlay::sequence<size_t> computeLeftEndpoints(int num_levels,
                                              size_t size_per_leaf,
                                              size_t remainder) {
  assert(num_levels < (int)split_points.size());

  parlay::sequence<size_t> to_sum(1 + split_points[num_levels].size());
  to_sum[0] = 0;
  parlay::parallel_for(0, split_points[num_levels].size(), [&](size_t i) {
    to_sum[i + 1] = size_per_leaf + ((split_points[num_levels][i] <= remainder) ? 1 : 0);
  });

  // to_sum -> holds 0 and then the number of points in each child
  //  - inclusive prefix sum - tells us endpoints for the children

  return parlay::scan_inclusive(to_sum);
}

// compute left indices of children at a given level with a given remainder
parlay::sequence<size_t> computeBottomChildSizes(int num_levels,
                                                 size_t size_per_leaf,
                                                 size_t remainder) {
  assert(num_levels < (int)split_points.size());

  parlay::sequence<size_t> to_sum(split_points[num_levels].size());
  parlay::parallel_for(0, split_points[num_levels].size(), [&](size_t i) {
    to_sum[i] =
        numNodesBottom(size_per_leaf + ((split_points[num_levels][i] <= remainder) ? 1 : 0));
  });

  // to_sum -> holds the number of nodes in each child
  //  - prefix sum - tells us where each child starts in memory

  return parlay::scan(to_sum).first;
}

std::vector<std::vector<size_t>> child_indices = {{}};
// child_indices[n]: the indices of the children of a tree with n levels, offset from the root
//  -> n levels => 2^{n-1} children, 2^n-1 total nodes
void initializeChildIndices(int num_levels) {
  assert(num_levels > 1);
  if (num_levels < (int)child_indices.size()) return;  // already done

  int level = (int)child_indices.size();  // the next level to start populating (because of the +1)
  child_indices.resize(num_levels + 1);   // allocate space for the levels

  // base case
  if (level == 1) {  // first time initializing
    child_indices[1] = {0};
    level++;
  }

  for (; level <= num_levels; level++) {
    // if(child_indices[level].size() > 0) continue; // skip completed
    auto bottom_num_levels = bottomNumLevels(level);
    auto top_num_levels = level - bottom_num_levels;

    auto top_tree_size = (1 << (top_num_levels)) - 1;
    auto bottom_tree_size = (1 << (bottom_num_levels)) - 1;

    auto num_bottom_trees = 1 << (top_num_levels);
    for (int subtree = 0; subtree < num_bottom_trees; subtree++) {
      for (const auto& child_idx : child_indices[bottom_num_levels]) {
        child_indices[level].push_back(top_tree_size + subtree * bottom_tree_size + child_idx);
      }
    }
  }
}

static bool initialized = false;
void initialize(__attribute__((unused)) int num_levels) {
  if (!initialized) {
    initializeChildIndices(24);
    initializeSplitPoints(24);
    initialized = true;
  }
}

} // End namespace batchKdTree