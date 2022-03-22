#ifndef LOGTREE_H
#define LOGTREE_H

#include "../cache-oblivious/cokdtree.h"
#include "../binary-heap-layout/bhlkdtree.h"
#include "../shared/macro.h"
#include "./buffer.h"

#ifdef PRINT_LOGTREE_TIMINGS
#include "common/get_time.h"
#endif

#ifdef LOGTREE_USE_BLOOM
#include "../shared/bloom.h"
#endif

/*#define DEBUG_MSG(str)             \
  do {                             \
    std::cout << str << std::endl; \
  } while (false)*/

#ifndef NDEBUG
template <int NUM_TREES>
static std::string treeMaskToString(int tree_mask) {
  std::stringstream ss;
  for (int i = NUM_TREES - 1; i >= 0; i--) {
    ss << ((tree_mask >> i) & 1);
  }
  return ss.str();
}
static std::string seq_to_str(const parlay::sequence<int>& s) {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < s.size(); i++) {
    if (i > 0) ss << ", ";
    ss << s[i];
  }
  ss << "]";
  return ss.str();
}
#endif

template <int NUM_TREES,         // the number of static trees
          int BUFFER_LOG2_SIZE,  // the size of the (dynamic) buffer tree
          int dim,
          class objT,
          bool parallel,  // defaults in kdtree.h forward decl
          bool coarsen>
class LogTree {
#if (LOGTREE_BUFFER == BHL_BUFFER)
  typedef BHL_KdTree<dim, objT, parallel, false> dynamicTree;
#else
  typedef LogTreeBuffer<dim, objT, parallel> dynamicTree;
#endif
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
  typedef CO_KdTree<dim, objT, parallel, coarsen> staticTree;
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
  typedef BHL_KdTree<dim, objT, parallel, coarsen> staticTree;
#endif
  typedef pargeo::point<dim> pointT;
#ifdef LOGTREE_USE_BLOOM
  typedef BloomFilter<dim> BloomFilterT;
#endif
  static const size_t BUFFER_SIZE = (1 << BUFFER_LOG2_SIZE);

  int tree_mask;  // represent whether the static trees are full or not
  dynamicTree buffer_tree;
#ifdef LOGTREE_USE_BLOOM
  BloomFilterT buffer_bloom_filter;
#endif

  // TODO: make sure this is on individual cache lines
  staticTree* static_trees;  // TODO: maybe use a vector?
#ifdef LOGTREE_USE_BLOOM
  BloomFilterT* static_bloom_filters;
#endif

  static inline int nth_tree_log2size(int n) {
    assert(n < NUM_TREES);
    return (n + BUFFER_LOG2_SIZE);
  }
  // TODO: make this size_t
  static inline size_t nth_tree_size(int n) { return 1 << nth_tree_log2size(n); }
  static inline bool nth_bit_set(int x, int n) { return (x >> n) % 2; }
  static inline void set_nth_bit(int& x, int n) { x |= (1 << n); }
  static inline void unset_nth_bit(int& x, int n) { x &= ~(1 << n); }

 public:
#ifdef ERASE_SEARCH_TIMES
  double total_search_time = 0;
  double total_bbox_time = 0;
  double total_leaf_time = 0;
#endif
  static constexpr bool coarsen_ = coarsen;
  LogTree()
      : tree_mask(0),
        buffer_tree(BUFFER_LOG2_SIZE, true)
#ifdef LOGTREE_USE_BLOOM
        ,
        buffer_bloom_filter(1 << BUFFER_LOG2_SIZE)
#endif
  {
    static_trees = (staticTree*)malloc(NUM_TREES * sizeof(staticTree));
#ifdef LOGTREE_USE_BLOOM
    static_bloom_filters = (BloomFilterT*)malloc(NUM_TREES * sizeof(BloomFilterT));
#endif
    for (int i = 0; i < NUM_TREES; i++) {
      auto curlog2size = nth_tree_log2size(i);
#if (PARTITION_TYPE == PARTITION_OBJECT_MEDIAN)
      new (&static_trees[i]) staticTree(curlog2size);
#elif (PARTITION_TYPE == PARTITION_SPATIAL_MEDIAN)
      new (&static_trees[i]) staticTree(curlog2size, false);
#endif
#ifdef LOGTREE_USE_BLOOM
      new (&static_bloom_filters[i]) BloomFilterT(1 << curlog2size);
#endif
    }
  }

  // hack to get treeTime to work
  LogTree(__attribute__((unused)) int x) : LogTree() {}
  template <class R>
  LogTree(const R& points) : LogTree() {
    insert(points);
  }

  ~LogTree() {
    for (int i = 0; i < NUM_TREES; i++) {
      static_trees[i].~staticTree();  // destroy placement-new objects

#ifdef LOGTREE_USE_BLOOM
      static_bloom_filters[i].~BloomFilterT();
#endif
    }
    free(static_trees);
#ifdef LOGTREE_USE_BLOOM
    free(static_bloom_filters);
#endif
  }

  // MODIFY -----------------------------------------
  // Memory management: We construct new point arrays in this function, and then move them into the
  // new trees.
  // TODO: maybe something relating to a single contiguous array will work
  // TODO: think about coarsening parallel base cases
  template <class R>  // TODO: use [parlay::Range] concept later
  void insert(const R& points) {
#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
    timer t("[Insert]");
#endif
    // compute number of moving elements in terms of buffers
    int full_buffers = (int)(points.size() / BUFFER_SIZE);
    int remainder = (int)(points.size() % BUFFER_SIZE);
    bool use_buffer = false;

    // check if buffer is involved
    if (remainder + buffer_tree.size() > BUFFER_SIZE) {
      full_buffers++;
      remainder = (remainder + buffer_tree.size()) - BUFFER_SIZE;
      use_buffer = true;
    }

    DEBUG_MSG("Inserting " << points.size() << " points");
    DEBUG_MSG("full buffers, remainder, use_buffer = " << full_buffers << ", " << remainder << ", "
                                                       << (use_buffer ? "true" : "false"));

    // simulate which trees to gather
    // <buffer, points_start, points_end, gather trees, new tree>
    parlay::sequence<std::tuple<bool, size_t, size_t, parlay::sequence<int>, int>> moves;
    auto cur_points_end = points.size();
    bool have_used_buffer = false;
    auto new_tree_mask = tree_mask + full_buffers;
    DEBUG_MSG(treeMaskToString<NUM_TREES>(tree_mask) << "\n"
                                                     << treeMaskToString<NUM_TREES>(new_tree_mask));

    for (int i = NUM_TREES - 1; i >= 0;) {
      if (nth_bit_set(new_tree_mask, i) && !nth_bit_set(tree_mask, i)) {
        DEBUG_MSG("New Tree: " << i);
        // tree [i] was not filled before but is now

        // gather all the smaller trees
        parlay::sequence<int> to_gather;
        size_t size_gathered = 0;
        int j = i - 1;
        while (j >= 0) {
          if (nth_bit_set(new_tree_mask, j)) {
            if (nth_bit_set(tree_mask, j)) {
              DEBUG_MSG(" - Skipping Tree (11)" << j);
            } else {
              break;
            }

          } else {
            if (nth_bit_set(tree_mask, j)) {
              DEBUG_MSG(" - Gathering Tree " << j);
              to_gather.push_back(j);
              size_gathered += nth_tree_size(j);
            } else {
              DEBUG_MSG(" - Skipping Tree (00)" << j);
            }
          }
          j--;
        }

        auto num_from_points = nth_tree_size(i) - size_gathered;
        full_buffers -= (num_from_points / BUFFER_SIZE);
        DEBUG_MSG("num_from_points: " << num_from_points);

        bool cur_uses_buffer = false;
        if (use_buffer && (full_buffers == 0)) {
          assert(!have_used_buffer);
          have_used_buffer = true;
          DEBUG_MSG(" - Gathering Buffer Tree ");
          cur_uses_buffer = true;
          num_from_points -= buffer_tree.size();  // we moved buffer tree in
        }

        auto cur_points_start = cur_points_end - num_from_points;
        DEBUG_MSG("MOVE: [uses_buffer, point start, point end, trees, new tree] = ["
                  << (cur_uses_buffer ? "true" : "false") << ", " << cur_points_start << ", "
                  << cur_points_end << ", " << seq_to_str(to_gather) << ", " << i << "]");
        moves.push_back(
            {cur_uses_buffer, cur_points_start, cur_points_end, std::move(to_gather), i});

        // reset counters
        cur_points_end = cur_points_start;
        i = j;
      } else {  // otherwise, should be the same
        assert(nth_bit_set(new_tree_mask, i) == nth_bit_set(tree_mask, i));
        i--;
      }
    }
    assert(cur_points_end == (size_t)remainder);

    // safety assertion, should remove later for release
    if (full_buffers != 0)
      throw std::runtime_error("Not enough trees! : full_buffers = " +
                               std::to_string(full_buffers));

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
    std::cout << "[Insert] Serial Computation: " << t.get_next() << "\n";
#endif
    // REBUILD THE TREES -----------------------------

    // need to serially empty buffer if it's used
    parlay::sequence<objT> buffer_points;
    if (have_used_buffer) {
      buffer_points.resize(buffer_tree.size());
      buffer_tree.moveElementsTo(buffer_points.cut(0, buffer_points.size()));

#ifdef LOGTREE_USE_BLOOM
      buffer_bloom_filter.clear();
#endif
    }

    // insert remainder into buffer
    auto fill_buff_f = [&]() {
#ifdef LOGTREE_USE_BLOOM
      parlay::par_do([&]() { buffer_tree.insert(points.cut(0, remainder)); },
                     [&]() { buffer_bloom_filter.insert(points.cut(0, remainder)); });
#else
      buffer_tree.insert(points.cut(0, remainder));
#endif
    };
    // use the simulated moves above to construct new trees
    auto rebuild_static_f = [&](size_t i) {
      const auto& [uses_buffer, points_start, points_end, trees, new_tree] = moves[i];
      auto num_points = points_end - points_start;

      // gather items
      parlay::sequence<objT> cur_items;

      // compute where each set of elements goes into [cur_items]
      parlay::sequence<size_t> gather_endpoints;
      gather_endpoints.resize(1 + trees.size() + 1 + (uses_buffer ? 1 : 0));
      int cur_idx = 0;
      gather_endpoints[cur_idx++] = 0;                          // left endpoint
      for (int j = 0; j < (int)trees.size(); j++, cur_idx++) {  // tree endpoints
        gather_endpoints[cur_idx] = static_trees[trees[j]].size() + gather_endpoints[cur_idx - 1];
      }
      gather_endpoints[cur_idx] = num_points + gather_endpoints[cur_idx - 1];
      cur_idx++;
      if (uses_buffer) {
        gather_endpoints[cur_idx] = buffer_points.size() + gather_endpoints[cur_idx - 1];
        cur_idx++;
      }
      cur_items.resize(gather_endpoints[cur_idx - 1]);  // the full size

      // construct the elements
      auto construct_points = [&](size_t idx) {
        if (idx == trees.size() + 1) {
          // move buffer
          assert(gather_endpoints[idx + 1] - gather_endpoints[idx] == buffer_points.size());
          auto move_buffer_f = [&](size_t j) {
            cur_items[j + gather_endpoints[idx]] = buffer_points[j];
          };

          if (parallel) {
            parlay::parallel_for(0, buffer_points.size(), move_buffer_f);
          } else {
            for (size_t j = 0; j < buffer_points.size(); j++)
              move_buffer_f(j);
          }
        } else if (idx == trees.size()) {
          // move points
          assert(gather_endpoints[idx + 1] - gather_endpoints[idx] == num_points);
          auto move_points_f = [&](size_t j) {
            cur_items[j + gather_endpoints[idx]] = points[j + points_start];
          };

          if (parallel) {
            parlay::parallel_for(0, num_points, move_points_f);
          } else {
            for (size_t j = 0; j < num_points; j++)
              move_points_f(j);
          }
        } else {
          // [0, num_trees) -> move a tree
          auto tree_idx = trees[idx];
          assert(((int)tree_idx < new_tree));
          assert(gather_endpoints[idx + 1] - gather_endpoints[idx] ==
                 static_trees[tree_idx].size());
          static_trees[tree_idx].moveElementsTo(
              cur_items.cut(gather_endpoints[idx], gather_endpoints[idx + 1]));
        }
      };

      if (parallel) {
        parlay::parallel_for(0, gather_endpoints.size() - 1, construct_points);
      } else {
        for (size_t idx = 0; idx < gather_endpoints.size() - 1; idx++)
          construct_points(idx);
      }
      //#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
      // if (new_tree == 16) {
      // std::cout << "(16) Point Gather: " << t.get_next() << "\n";
      //}
      //#endif

      // construct the new tree
      assert(static_trees[new_tree].empty());
      DEBUG_MSG("CONSTRUCTING TREE[" << new_tree << "]: " << cur_items.size() << " items");

#ifdef LOGTREE_USE_BLOOM
#ifdef BLOOM_FILTER_BUILD_COPY
      auto items_copy = cur_items;
      parlay::par_do([&]() { static_trees[new_tree].build(std::move(cur_items)); },
                     [&]() { static_bloom_filters[new_tree].build(items_copy); });
#else
      static_trees[new_tree].build(std::move(cur_items));
      static_bloom_filters[new_tree].build(static_trees[new_tree].items);
#endif
#else
      static_trees[new_tree].build(std::move(cur_items));
#endif

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
      std::cout << "[Insert] Tree[" << new_tree << "] Construction Time: " << t.get_next() << "\n";
#endif
      //#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
      // if (new_tree == 16) {
      // std::cout << "(16) Tree Build: " << t.get_next() << "\n";
      //}
      //#endif
    };

    if (parallel) {
      parlay::parallel_for(
          0,
          moves.size() + 1,
          [&](size_t i) {
            if (i == moves.size())
              fill_buff_f();
            else
              rebuild_static_f(i);
          },
          1);
    } else {
      fill_buff_f();
      // rebuild the static trees
      for (int i = 0; i < (int)moves.size(); i++) {
        rebuild_static_f(i);
      }
    }

    // update tree mask
    tree_mask = new_tree_mask;

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_INSERT_TIMINGS)
    std::cout << "[Insert] Build: " << t.get_next() << "\n";
    t.reportTotal("Total");
    t.stop();
#endif
  }

  template <bool bulk, class R>
  void erase(const R& points) {
    constexpr int BUFFER_TREE_IDX = -1;
    auto tree_ids = gatherFullTrees();
    // PHASE 1: Erase points from trees ----------------------
    auto erase_from_tree = [&](size_t tid) {
      // DEBUG_MSG("Erasing from tree " << i << "/" << NUM_TREES);
      auto i = tree_ids[tid];

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_DELETE_TIMINGS)
      timer t("[Delete]");
#endif

      parlay::sequence<pointT> to_erase;
#ifdef LOGTREE_USE_BLOOM
      if (i == BUFFER_TREE_IDX)
        to_erase = buffer_bloom_filter.filter(points);
      else
        to_erase = static_bloom_filters[i].filter(points);
#else
      to_erase = points;
#endif

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_DELETE_TIMINGS)
      std::stringstream ss;
      ss << "[Delete] Filter[" << i << "]: " << t.get_next() << "\n";
      ss << "         -> " << to_erase.size() << " / " << points.size() << std::endl;
#endif

      if (bulk) {
        if (i == BUFFER_TREE_IDX) {
          buffer_tree.template bulk_erase<false>(to_erase);
#ifdef ERASE_SEARCH_TIMES
          total_search_time += buffer_tree.total_search_time;
#endif
        } else {
          static_trees[i].template bulk_erase<false>(to_erase);
#ifdef ERASE_SEARCH_TIMES
          total_search_time += static_trees[i].total_search_time;
#endif
        }
      } else {
        if (i == BUFFER_TREE_IDX) {
          buffer_tree.template erase<true>(to_erase);
        } else {
          static_trees[i].template erase<true>(to_erase);
        }
      }

#if defined(PRINT_LOGTREE_TIMINGS) && defined(PRINT_DELETE_TIMINGS)
      ss << "[Delete] Delete[" << i << "]: " << t.get_next() << "\n";
      printf("%s", ss.str().c_str());
#endif
    };

    if (parallel) {
      parlay::parallel_for(0, tree_ids.size(), erase_from_tree, 1);
    } else {
      for (size_t i = 0; i < tree_ids.size(); i++) {
        erase_from_tree(i);
      }
    }

    // PHASE 2: Collect all depleted trees
    // compute depleted trees
    parlay::sequence<size_t> gather_points;
    parlay::sequence<int> depleted_trees;
    auto new_tree_mask = tree_mask;
    gather_points.push_back(0);  // initialize
    for (int i = 0; i < NUM_TREES; i++) {
      if (static_trees[i].size() <= nth_tree_size(i) / 2) {
        // need to push down
        depleted_trees.push_back(i);
        gather_points.push_back(gather_points.back() + static_trees[i].size());
        unset_nth_bit(new_tree_mask, i);
      }
    }
    tree_mask = new_tree_mask;

    // gather depleted trees
    parlay::sequence<objT> points_to_move(gather_points.back());
    auto gather_tree = [&](size_t i) {
      auto tree_idx = depleted_trees[i];
      assert(static_trees[tree_idx].size() == gather_points[i + 1] - gather_points[i]);
      static_trees[tree_idx].moveElementsTo(
          points_to_move.cut(gather_points[i], gather_points[i + 1]));
    };
    if (parallel) {
      parlay::parallel_for(0, depleted_trees.size(), gather_tree);
    } else {
      for (size_t i = 0; i < depleted_trees.size(); i++)
        gather_tree(i);
    }

    // reinsert them
    insert(points_to_move);
  }

  template <class R>
  void bulk_erase(const R& points) {
    erase<true, R>(points);
  }

  /*
  // TODO: deduplicate points in bulk erase
  template <bool bulk>
  void erase_old(const parlay::sequence<objT>& points) {
    // PHASE 1: Erase points from trees ----------------------
    parlay::sequence<int> delete_counts(NUM_TREES + 1);

    auto erase_from_tree = [&](size_t i) {
      auto orig_size = (i == NUM_TREES) ? buffer_tree.size() : static_trees[i].size();
      if (bulk) {
        if (i == NUM_TREES) {
          buffer_tree.template bulk_erase<false>(points);
        } else {
          static_trees[i].template bulk_erase<false>(points);
        }
      } else {
        if (i == NUM_TREES) {
          buffer_tree.template erase<true>(points);
        } else {
          static_trees[i].template erase<true>(points);
        }
      }
      auto new_size = (i == NUM_TREES) ? buffer_tree.size() : static_trees[i].size();
      delete_counts[i] = orig_size - new_size;
    };

    if (parallel) {
      parlay::parallel_for(0, NUM_TREES + 1, erase_from_tree, 1);
    } else {
      for (int i = 0; i < NUM_TREES + 1; i++) {
        erase_from_tree(i);
      }
    }

    // PHASE 2: Compute tree movements -----------------------
    // simulate + figure out movements
    static const int NORMAL_ERASE = -3;
    static const int MOVE_DOWN = -2;
    static const int DYNAMIC_BUFFER = -1;

    int move_from_buffer = 0;
    parlay::sequence<std::pair<int, int>> changes;

    // dynamic tree delete
#ifndef NDEBUG
    if (delete_counts[NUM_TREES] > 0) changes.emplace_back(NORMAL_ERASE, DYNAMIC_BUFFER);
#endif

    // TODO: switch this to only touch nonempty trees; consider unrolling tree 0 out
    // TODO: speed up the case where a tree is completely emptied
    auto new_tree_mask = tree_mask;
    for (int i = 0; i < NUM_TREES; i++) {
      if (!nth_bit_set(new_tree_mask, i)) continue;
      if (delete_counts[i] == 0) continue;

      auto new_size = static_trees[i].size();
      if (new_size > nth_tree_size(i) / 2) {
#ifndef NDEBUG
        changes.emplace_back(NORMAL_ERASE, i);  // only record a normal erase for debugging
#endif
      } else {
        if (i == 0) {                               // special case: have to think about buffer tree
          auto cutoff = BUFFER_SIZE - new_size;     // how much tree 0 needs to be full
          if ((int)buffer_tree.size() <= cutoff) {  // buffer doesn't have enough -> move 0 down
            changes.emplace_back(MOVE_DOWN, i);
            unset_nth_bit(new_tree_mask, i);
          } else {  // buffer has enough -> move enough up to fill 0
            changes.emplace_back(DYNAMIC_BUFFER, 0);
            assert(move_from_buffer == 0);
            move_from_buffer = cutoff;
          }
        } else {
          if (!nth_bit_set(new_tree_mask, i - 1)) {
            changes.emplace_back(MOVE_DOWN, i);  // just push it down
            unset_nth_bit(new_tree_mask, i);     // move down to i-1
            set_nth_bit(new_tree_mask, i - 1);
          } else {
            if (changes.size() == 0 || changes.back().first <= MOVE_DOWN ||
                changes.back().second != i - 1) {
              changes.emplace_back(i - 1, i);
            } else {  // need to cascade
              changes.back().second = i;
            }
            unset_nth_bit(new_tree_mask, i - 1);  // move i-1 up to i
          }
        }
      }
    }

#ifndef NDEBUG
    DEBUG_MSG("DELETE: " << tree_mask << " -> " << new_tree_mask);
    for (const auto& p : changes) {
      if (p.first == NORMAL_ERASE) {
        int num_points;
        if (p.second == DYNAMIC_BUFFER) {
          num_points = delete_counts[NUM_TREES];
        } else
          num_points = delete_counts[p.second];
        std::cout << "ERASE[" << p.second << "]: " << num_points << "/"
                  << ((p.second >= 0) ? nth_tree_size(p.second) : BUFFER_SIZE) << " points"
                  << std::endl;
      } else if (p.first == MOVE_DOWN) {
        auto new_size = static_trees[p.second].size();
        auto cur_size = new_size + delete_counts[p.second];
        auto full_size = nth_tree_size(p.second);
        std::cout << "MOVE_DOWN[" << p.second << " -> " << p.second - 1 << "]: " << cur_size
                  << " -> " << new_size << " / " << full_size << " points" << std::endl;
      } else {
        std::cout << "MOVE_UP[" << p.first << " -> " << p.second << "]" << std::endl;
      }
    }
#endif

    // PHASE 3: Carry out simulated changes -----------------------------

    auto perform_change = [&](size_t i) {
      const auto& p = changes[i];
      if (p.first == NORMAL_ERASE) {
        // do nothing
      } else if (p.first == MOVE_DOWN) {
        assert((p.second >= 0) && (p.second < NUM_TREES));
        // pull the elements out
        parlay::sequence<objT> items_to_move(static_trees[p.second].size());
        [[maybe_unused]] auto num_moved =
            static_trees[p.second].moveElementsTo(items_to_move.cut(0, items_to_move.size()));
        assert(num_moved == items_to_move.size());
        if (p.second == 0) {
          assert(buffer_tree.size() + items_to_move.size() <= BUFFER_SIZE);
          const auto& const_items = items_to_move;
          buffer_tree.insert(const_items.cut(0, const_items.size()));
        } else {
          assert(static_trees[p.second - 1].empty());
          static_trees[p.second - 1].build(std::move(items_to_move));
        }
      } else {  // need to gather all the points
        // PHASE 1: compute size to move (serial because <= NUM_TREES trees total) -------
        parlay::sequence<size_t> move_offsets(p.second - p.first + 2);
        move_offsets[0] = 0;

        for (auto idx = p.first; idx <= p.second; idx++) {
          size_t to_add;
          if (idx == DYNAMIC_BUFFER) {
            to_add = move_from_buffer;
          } else {
            to_add = static_trees[idx].size();
          }
          move_offsets[idx - p.first + 1] = move_offsets[idx - p.first] + to_add;
        }

        // allocate space
        parlay::sequence<objT> items_to_move(move_offsets.back());

        // PHASE 2: gather all the elements ----------------
        auto gather_elements = [&](size_t i) {
          int idx = (int)i - 1;
          if (idx == DYNAMIC_BUFFER) {
            // get buffer elements
            parlay::sequence<objT> buffer_items(buffer_tree.size());
            buffer_tree.moveElementsTo(buffer_items.cut(0, buffer_items.size()));

            auto leave_behind = [&]() {
              // leave behind the extra elements
              const auto& const_buffer_items = buffer_items;
              buffer_tree.insert(const_buffer_items.cut(move_from_buffer, buffer_items.size()));
            };

            auto keep = [&]() {
              // pick up the pieces to move
              // TODO: could probably take more, depending on sizes of later trees
              parlay::parallel_for(
                  0, move_from_buffer, [&](size_t i) { items_to_move[i] = buffer_items[i]; });
            };

            // redistribute the buffer tree
            if (parallel) {
              parlay::par_do(leave_behind, keep);
            } else {
              leave_behind();
              keep();
            }

          } else {
            // gather the points
            assert(move_offsets[idx - p.first + 1] ==
                   move_offsets[idx - p.first] + static_trees[idx].size());
            static_trees[idx].moveElementsTo(
                items_to_move.cut(move_offsets[idx - p.first], move_offsets[idx - p.first + 1]));
          }
        };

        if (parallel) {
          parlay::parallel_for(p.first + 1, p.second + 2, gather_elements);
        } else {
          for (size_t i = (size_t)(p.first + 1); i < (size_t)p.second + 2; i++) {
            gather_elements(i);
          }
        }

        // move the elements
        static_trees[p.second].build(std::move(items_to_move));
      }
    };

    if (parallel) {
      parlay::parallel_for(0, changes.size(), perform_change);
    } else {
      for (size_t i = 0; i < changes.size(); i++)
        perform_change(i);
    }

    // update
    tree_mask = new_tree_mask;
  }
  */

  // QUERY -----------------------------------------
  bool contains(const objT& p) const {
    if (parallel) {
      parlay::sequence<bool> res(NUM_TREES + 1);
      parlay::parallel_for(0, NUM_TREES + 1, [&](size_t i) {
        if (i == NUM_TREES) {
          res[i] = buffer_tree.contains(p);
        } else {
          res[i] = static_trees[i].contains(p);
        }
      });

      // extract the result
      for (const auto& b : res)
        if (b) return true;
      return false;
    } else {
      if (buffer_tree.contains(p)) return true;
      for (int i = 0; i < NUM_TREES; i++)
        if (static_trees[i].contains(p)) return true;
      return false;
    }
  }

  parlay::sequence<objT> orthogonalQuery(const objT& qMin, const objT& qMax) const {
    if (parallel) {
      parlay::sequence<parlay::sequence<objT>> res(NUM_TREES + 1);
      parlay::parallel_for(0, NUM_TREES + 1, [&](size_t i) {
        if (i == NUM_TREES) {
          res[i] = buffer_tree.orthogonalQuery(qMin, qMax);
        } else {
          res[i] = static_trees[i].orthogonalQuery(qMin, qMax);
        }
      });

      // result size
      std::vector<size_t> offsets(NUM_TREES + 2);
      offsets[0] = 0;
      for (int i = 1; i < NUM_TREES + 2; i++) {
        offsets[i] = offsets[i - 1] + res[i - 1].size();
      }

      // reduce the result
      parlay::sequence<objT> ret(offsets[NUM_TREES + 1]);
      parlay::parallel_for(0, res.size(), [&](size_t i) {
        assert(offsets[i + 1] - offsets[i] == res[i].size());
        parlay::parallel_for(0, res[i].size(), [&](size_t j) { ret[j + offsets[i]] = res[i][j]; });
      });
      return ret;
    } else {
      parlay::sequence<objT> ret;
      auto r = buffer_tree.orthogonalQuery(qMin, qMax);
      ret.insert(ret.begin(), r.begin(), r.end());
      for (int i = 0; i < NUM_TREES; i++) {
        r = static_trees[i].orthogonalQuery(qMin, qMax);
        ret.insert(ret.begin() + ret.size(), r.begin(), r.end());
      }
      return ret;
    }
  }

  parlay::sequence<int> gatherFullTrees() const {
    constexpr int BUFFER_TREE_IDX = -1;
    // gather full trees
    parlay::sequence<int> tree_ids;
    if (!buffer_tree.empty()) {
      tree_ids.push_back(BUFFER_TREE_IDX);
    }
    for (int i = 0; i < NUM_TREES; i++) {
      if (nth_bit_set(tree_mask, i)) tree_ids.push_back(i);
    }
    return tree_ids;
  }

  template <bool update = false, bool recurse_sibling = false>
  parlay::sequence<const pointT*> knn3(const parlay::sequence<objT>& queries, int k) const {
#ifdef PRINT_LOGTREE_TIMINGS
    timer t;
#endif
    constexpr int BUFFER_TREE_IDX = -1;
    auto tree_ids = gatherFullTrees();

    // result buffer
    parlay::sequence<const pointT*> res(k * queries.size());
    auto res_slice = res.head(res.size());

    // knn buffer
    auto out_size = (2 * k * queries.size());
    parlay::sequence<knnBuf::elem<const pointT*>> out(out_size);
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN3] Space Allocation: " << t.get_next() << "\n";
#endif

    // call knn on one tree at a time, but parallel within tree
    auto run_on_tree = [&](size_t i, auto out_slice, bool preload) {
#ifdef PRINT_LOGTREE_TIMINGS
      timer t1;
#endif
      // get the id of the points
      auto tree_id = tree_ids[i];

      // call knn on this tree
      if (tree_id == BUFFER_TREE_IDX) {
        buffer_tree.template knn<false, update, recurse_sibling>(
            queries, out_slice, res_slice, k, preload);
      } else {
        static_trees[tree_id].template knn<false, update, recurse_sibling>(
            queries, out_slice, res_slice, k, preload);
      }
#ifdef PRINT_LOGTREE_TIMINGS
      std::cout << "[KNN3] Tree " << tree_id << " Query Time: " << t1.get_next() << "\n";
#endif
    };

    auto out_slice = out.cut(0, out.size());  // use the same slice at every iteration
    for (int i = (int)tree_ids.size() - 1; i >= 0; i--)
      run_on_tree(i, out_slice, i != ((int)tree_ids.size() - 1));
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN3] Query Time: " << t.get_next() << "\n";
#endif

    // combine results
    for (size_t i = 0; i < queries.size(); i++) {
      for (int g = 0; g < k; g++) {
        res[i * k + g] = out[i * 2 * k + g].entry;
      }
    }
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN3] End: " << t.get_next() << "\n";
    t.reportTotal("[KNN3] Total");
    t.stop();
#endif

    return res;
  }

  template <bool update = false, bool recurse_sibling = false>
  parlay::sequence<const pointT*> knn2(const parlay::sequence<objT>& queries, int k) const {
    constexpr int BUFFER_TREE_IDX = -1;
    auto tree_ids = gatherFullTrees();

    // result buffer
    parlay::sequence<const pointT*> res(k * queries.size());
    auto res_slice = res.head(res.size());

    // knn buffer
    auto out_size = (2 * k * queries.size());
    parlay::sequence<knnBuf::elem<const pointT*>> out(out_size);
    auto out_slice = out.head(out.size());

    auto run_on_point = [&](size_t i) {
      // knn point i through all the trees
      for (size_t j = 0; j < tree_ids.size(); j++) {
        auto tree_id = tree_ids[j];
        auto preload = j > 0;  // buffer is full after first tree
        if (tree_id == BUFFER_TREE_IDX) {
          buffer_tree.template knnSinglePoint<false, update, recurse_sibling>(
              queries[i], i, out_slice, res_slice, k, preload);
        } else {
          static_trees[tree_id].template knnSinglePoint<false, update, recurse_sibling>(
              queries[i], i, out_slice, res_slice, k, preload);
        }
      }
      // gather results
      for (int j = 0; j < k; j++) {
        res[i * k + j] = out_slice[(i * 2 * k) + j].entry;
      }
    };

    if (parallel) {
      parlay::parallel_for(0, queries.size(), run_on_point);
    } else {
      for (size_t i = 0; i < queries.size(); i++) {
        run_on_point(i);
      }
    }

    return res;
  }

  template <bool update = false, bool recurse_sibling = false>
  parlay::sequence<const pointT*> knn(const parlay::sequence<objT>& queries, int k) const {
#ifdef PRINT_LOGTREE_TIMINGS
    timer t;
#endif
    constexpr int BUFFER_TREE_IDX = -1;
    auto tree_ids = gatherFullTrees();

    // result buffer
    parlay::sequence<const pointT*> res(k * queries.size());
    auto res_slice = res.head(res.size());

    // knn buffer
    auto out_size = (2 * k * queries.size());
    parlay::sequence<knnBuf::elem<const pointT*>> out(out_size * (parallel ? tree_ids.size() : 1));
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN] Space Allocation: " << t.get_next() << "\n";
#endif

    // call knn in parallel on all trees
    auto run_on_tree = [&](size_t i, auto out_slice, bool preload) {
#ifdef PRINT_LOGTREE_TIMINGS
      timer t1;
#endif
      // get the id of the points
      auto tree_id = tree_ids[i];

      // call knn on this tree
      if (tree_id == BUFFER_TREE_IDX) {
        buffer_tree.template knn<false, update, recurse_sibling>(
            queries, out_slice, res_slice, k, preload);
      } else {
        static_trees[tree_id].template knn<false, update, recurse_sibling>(
            queries, out_slice, res_slice, k, preload);
      }
#ifdef PRINT_LOGTREE_TIMINGS
      std::cout << "[KNN] Tree " << tree_id << " Query Time: " << t1.get_next() << "\n";
#endif
    };

    if (parallel) {
      parlay::parallel_for(0, tree_ids.size(), [&](size_t i) {
        auto out_slice = out.cut(i * out_size, (i + 1) * out_size);
        run_on_tree(i, out_slice, false);
      });
    } else {
      auto out_slice = out.cut(0, out.size());  // use the same slice at every iteration
      for (size_t i = 0; i < tree_ids.size(); i++)
        run_on_tree(i, out_slice, i > 0);
    }
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN] Query Time: " << t.get_next() << "\n";
#endif

    // combine results
    if (parallel) {
      parlay::parallel_for(0, queries.size(), [&](size_t i) {
        auto buf = knnBuf::buffer<const pointT*>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
        buf.ptr = k;
        for (size_t j = 1; j < tree_ids.size(); j++) {
          auto start_elem = out.begin() + (j * out_size + i * 2 * k);
          for (int g = 0; g < k; g++) {
            buf.insert(*(start_elem + g));
          }
        }
        for (int g = 0; g < k; g++) {
          res[i * k + g] = buf[g].entry;
        }
      });
    } else {
      for (size_t i = 0; i < queries.size(); i++) {
        for (int g = 0; g < k; g++) {
          res[i * k + g] = out[i * 2 * k + g].entry;
        }
      }
    }
#ifdef PRINT_LOGTREE_TIMINGS
    std::cout << "[KNN] End: " << t.get_next() << "\n";
    t.reportTotal("[KNN] Total");
    t.stop();
#endif

    return res;
  }

  parlay::sequence<const pointT*> dualKnnBase(const KdTree<dim, objT, parallel, coarsen>& queryTree,
                                              int k) const {
#ifdef PRINT_LOGTREE_TIMINGS
    timer t;
#endif

    constexpr int BUFFER_TREE_IDX = -1;
    auto tree_ids = gatherFullTrees();

    // result buffer
    parlay::sequence<const pointT*> res(k * queryTree.size());

    // knn buffer
    auto out_size = (2 * k * queryTree.size());
#if (DUAL_KNN_MODE == DKNN_NONATOMIC_LEAF)
    auto multiplier = 1;  // in non-atomic case, run on one tree at a time -> only one buffer
#else
    auto multiplier =
        (parallel ? tree_ids.size() : 1);  // in atomic/array case, have separate buffers for each
#endif
    parlay::sequence<knnBuf::elem<const pointT*>> out(out_size * multiplier);
    parlay::sequence<knnBuf::buffer<const pointT*>> bufs(queryTree.size() * multiplier);

#if (DUAL_KNN_MODE == DKNN_NONATOMIC_LEAF)
    // in non-atomic case, have to initialize buffers ahead of time
    //  -> in other cases, initialize buffers in each tree's thread
    auto setup_buffer = [&](size_t i) {
      bufs[i] = knnBuf::buffer<const pointT*>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
    };
    if (parallel) {
      parlay::parallel_for(0, bufs.size(), setup_buffer);
    } else {
      for (size_t i = 0; i < bufs.size(); i++) {
        setup_buffer(i);
      }
    }
#endif

    // call knn in parallel on all trees
    auto run_on_tree = [&](size_t i, auto& buf_slice) {
#if (DUAL_KNN_MODE != DKNN_NONATOMIC_LEAF)
      // set up the buffers
      auto setup_buffer = [&](size_t j) {
        auto idx = i * queryTree.size() + j;
        bufs[idx] = knnBuf::buffer<const pointT*>(k, out.cut(idx * 2 * k, (idx + 1) * 2 * k));
      };
      if (parallel) {
        parlay::parallel_for(0, queryTree.size(), setup_buffer);
      } else if (i == 0) {  // only setup once in serial case
        for (size_t j = 0; j < queryTree.size(); j++) {
          setup_buffer(j);
        }
      }
#endif
      // get the id of the points
      auto tree_id = tree_ids[i];

      // call knn on this tree

#if (DUAL_KNN_MODE == DKNN_ARRAY)
      parlay::sequence<double> dualKnnDists(queryTree.num_nodes(),
                                            std::numeric_limits<double>::max());
#endif
      if (tree_id == BUFFER_TREE_IDX) {
        //#if (LOGTREE_BUFFER == BHL_BUFFER)
        //#if (DUAL_KNN_MODE == DKNN_ARRAY)
        // DualKnnHelper(queryTree.unsafe_root(),
        // buffer_tree.root(),
        // queryTree,
        // dualKnnDists,
        // buffer_tree,
        // buf_slice);
        //#else
        // DualKnnHelper(
        // queryTree.unsafe_root(), buffer_tree.root(), queryTree, buffer_tree, buf_slice);
        //#endif
        //#else
        //timer t; // TODO
#if (LOGTREE_BUFFER == ARR_BUFFER)
        buffer_tree.knn(queryTree.items, buf_slice);
#else
        buffer_tree.template knn<false, false>(queryTree.items, buf_slice);
#endif
#ifdef PRINT_LOGTREE_TIMINGS
        std::cout << "[DKNN] SMALL KNN: " << t.get_next() << "\n";
#endif
#if (DUAL_KNN_MODE == DKNN_ARRAY)
        queryTree.nodes[0].updateDualDist(
            buf_slice, queryTree.items.begin(), queryTree.nodes, dualKnnDists);
#else
        queryTree.nodes[0].updateDualDist(buf_slice, queryTree.items.begin());
#endif
#ifdef PRINT_LOGTREE_TIMINGS
        std::cout << "[DKNN] QTREE UPDATE: " << t.get_next() << "\n";
#endif
        //#endif
      } else {
#if (DUAL_KNN_MODE == DKNN_ARRAY)
        DualKnnHelper(queryTree.unsafe_root(),
                      static_trees[tree_id].root(),
                      queryTree,
                      dualKnnDists,
                      static_trees[tree_id],
                      buf_slice);
#else
        DualKnnHelper(queryTree.unsafe_root(),
                      static_trees[tree_id].root(),
                      queryTree,
                      static_trees[tree_id],
                      buf_slice);
#endif
      }
    };

#if (DUAL_KNN_MODE != DKNN_NONATOMIC_LEAF)
    if (parallel) {
      parlay::parallel_for(0, tree_ids.size(), [&](size_t i) {
        auto buf_slice = bufs.cut(i * queryTree.size(), (i + 1) * queryTree.size());
        run_on_tree(i, buf_slice);
      });
    } else {
#endif
      auto buf_slice = bufs.cut(0, bufs.size());  // use the same slice at every iteration
#ifdef PRINT_LOGTREE_TIMINGS
      std::cout << "[DKNN] Setup: " << t.get_next() << std::endl;
#endif
      for (size_t i = 0; i < tree_ids.size(); i++) {
        run_on_tree(i, buf_slice);
#ifdef PRINT_LOGTREE_TIMINGS
        auto nth_tree_cur_size = [&](int id) {
          if (id == BUFFER_TREE_IDX) {
            return buffer_tree.size();
          } else {
            assert(id < NUM_TREES);
            assert(id >= 0);
            return static_trees[id].size();
          }
        };
        auto tree_id = tree_ids[i];
        std::cout << "[DKNN] Tree #" << tree_id << "(" << nth_tree_cur_size(tree_id)
                  << " pts): " << t.get_next() << std::endl;
#endif
      }
#if (DUAL_KNN_MODE != DKNN_NONATOMIC_LEAF)
    }
#endif

    // combine results
    if (parallel) {
#if (DUAL_KNN_MODE != DKNN_NONATOMIC_LEAF)
      parlay::parallel_for(0, queryTree.size(), [&](size_t i) {
        auto& buf = bufs[i];
        for (size_t j = 1; j < tree_ids.size(); j++) {
          auto start_elem = out.begin() + (j * out_size + i * 2 * k);
          for (int g = 0; g < k; g++) {
            buf.insert(*(start_elem + g));
          }
        }
        buf.keepK();
        for (int g = 0; g < k; g++) {
          res[i * k + g] = buf[g].entry;
        }
      });
#else
      parlay::parallel_for(0, queryTree.size(), [&](size_t i) {
        for (int g = 0; g < k; g++) {
          res[i * k + g] = bufs[i][g].entry;
        }
      });
#endif
    } else {
      for (size_t i = 0; i < queryTree.size(); i++) {
        for (int g = 0; g < k; g++) {
          res[i * k + g] = bufs[i][g].entry;
        }
      }
    }

    return res;
  }

  // DEBUG
  int getTreeMask() const { return tree_mask; }

  // TODO: can make this better by tracking as inserts/deletes are done
  size_t size() const {
    size_t res = buffer_tree.size();
    for (int i = 0; i < NUM_TREES; i++) {
      if (nth_bit_set(tree_mask, i)) res += static_trees[i].size();
    }
    return res;
  }

  void print(int tree_idx) const {
    if (tree_idx < 0 || tree_idx >= NUM_TREES) throw std::runtime_error("tree_idx out of bounds!");
    static_trees[tree_idx].print();
  }
};

#endif  // LOGTREE_H
