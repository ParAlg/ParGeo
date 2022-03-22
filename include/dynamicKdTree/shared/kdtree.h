#ifndef KDTREE_H
#define KDTREE_H

#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "pargeo/point.h"

#include "kdnode.h"
#include "utils.h"
#include "knnbuffer.h"
#include "box.h"
#include "macro.h"

#ifdef ALL_USE_BLOOM
#include "bloom.h"
#endif

#ifdef PRINT_KDTREE_TIMINGS
#include "common/get_time.h"
#endif

// TODO: refactor so that there's static methods for checking if this is empty, constructing empty
// struct
// TODO: refactor this into an opaque substruct of KdTree, so that it can't be inspected outside
// this file.
// TODO: refactor this to use pointers, rather than idx
struct FoundPoint {
  static const long int NOT_FOUND = -1;
  long int idx;
  long int parent_idx;
  long int gparent_idx;
  int point_idx;  // index of point in [subtree_items]
};
[[maybe_unused]] static std::ostream &operator<<(std::ostream &os, FoundPoint &fp) {
  os << "{node_idx=" << fp.idx << ", parent_idx=" << fp.parent_idx
     << ", gparent_idx=" << fp.gparent_idx << ", point_idx=" << fp.point_idx << "}";
  return os;
}

// forward declare for dualknn
template <int NUM_TREES,
          int BUFFER_LOG2_SIZE,
          int dim,
          class objT,
          bool parallel = false,
          bool coarsen = false>
class LogTree;

template <int dim, class objT, bool parallel = false, bool coarsen = false>
class KdTree {
  template <int _NUM_TREES,
            int _BUFFER_LOG2_SIZE,
            int _dim,
            class _objT,
            bool _parallel,
            bool _coarsen>
  friend class LogTree;

 protected:
  typedef int intT;
  typedef double floatT;
  typedef pargeo::point<dim> pointT;
  typedef kdNode<dim, objT, parallel> nodeT;

#ifdef PRINT_KDTREE_TIMINGS
  inline void reset_timer() {
    timer_.reset();
    timer_.start();
  }

  inline void mark_time(const std::string &s) { timer_.report(timer_.get_next(), s); }
  inline void total_time() { timer_.reportTotal("Total"); }
#endif

  // TODO: this is probably not great for caching; also, it doesn't get zeroed out in release mode
  // => probably just want to get rid of it
  // nodeT **parents;
  nodeT *nodes;

  size_t cur_size;        // current number of nodes
  size_t build_size;      // the number of nodes it was built with
  const size_t max_size;  // the maximum size for this tree

  parlay::sequence<bool> present;
  parlay::sequence<objT> items;

#ifdef PRINT_KDTREE_TIMINGS
  timer timer_;
#endif

#ifdef ALL_USE_BLOOM
  typedef BloomFilter<dim> BloomFilterT;
  BloomFilterT bloom_filter;
#endif

 public:
#ifdef ERASE_SEARCH_TIMES
  double total_search_time;
  double total_bbox_time;
  double total_leaf_time;
#endif
  KdTree() = delete;
  KdTree(int log2size)
      : max_size(1UL << log2size)
#ifdef PRINT_KDTREE_TIMINGS
        ,
        timer_("KdTree")
#endif
#ifdef ALL_USE_BLOOM
        ,
        bloom_filter(1UL << log2size)
#endif
  {
#ifdef ERASE_SEARCH_TIMES
    total_search_time = 0;
    total_bbox_time = 0;
    total_leaf_time = 0;
#endif
    present = parlay::sequence<bool>(max_size);

    // TODO: use new[] for type safety
    nodes = (nodeT *)malloc((2 * max_size - 1) * sizeof(nodeT));

    // parents = (nodeT **)malloc((2 * max_size - 1) * sizeof(nodeT *));
    // parents[0] = nullptr;  // root

    clear();
  }

  // for convenience; construct the next largest tree
  template <class R>
  KdTree(const R &points) : KdTree((int)std::ceil(std::log2(points.size()))) {}

  ~KdTree() { free(nodes); }

  // MODIFY -----------------------------------------
  /*!
   * Clear out the contents of this tree
   */
  void clear() {
    cur_size = 0;
    build_size = 0;
    if (parallel) {
      // TODO: make this [present] reset better?
      // TODO: this is probably wrong - should reset to false!
      parlay::parallel_for(0, present.size(), [&](size_t i) { present[i] = true; });
    } else {
      for (size_t i = 0; i < present.size(); i++) {
        present[i] = true;
      }
    }
#ifdef ALL_USE_BLOOM
    bloom_filter.clear();
#endif
#ifndef NDEBUG
    // mark all the nodes as empty again, only for debugging purposes
    parlay::parallel_for(0, 2 * max_size - 1, [&](size_t i) { nodes[i].setEmpty(); });
    // parlay::parallel_for(0, 2 * n - 1, [&](size_t i) { parents[i] = nullptr; });
#endif
  }

  /*!
   * Move the elements of the tree and pack them into [dest]. Clear the tree.
   */
  size_t moveElementsTo(parlay::slice<objT *, objT *> dest) {
    size_t ret;
    assert(dest.size() >= size());
    if (parallel) {
      ret = parlay::pack_into(parlay::slice(items.begin(), items.begin() + build_size),
                              parlay::slice(present.begin(), present.begin() + build_size),
                              dest);
    } else {
      ret = 0;
      for (size_t i = 0; i < build_size; i++) {
        if (present[i]) {
          dest[ret++] = items[i];
        }
      }
    }
    clear();
    return ret;
  }

  // QUERY --------------------------------------------
  // return (s,e) where the node represents items [s,e) in the underlying [items] array
  std::pair<int, int> getNodeValueIdx(const nodeT *n) const {
    return {n->getStartValue() - items.begin(), n->getEndValue() - items.begin()};
  }
  // int getLeafValueIdx(const nodeT *n) const {
  // assert(n->isLeaf());
  // return getNodeValueIdx(n).first;
  //}

  bool contains(objT &p) const {
    if (empty()) return false;

    auto node = nodes;
    while (!node->isLeaf()) {  // walk down to leaf
      if (p.at(node->getSplitDimension()) < node->getSplitValue())
        node = node->getLeft();
      else
        node = node->getRight();
    }

    assert(node->isLeaf());
    // check in the leaf for the point
    for (auto it = node->getStartValue(); it != node->getEndValue(); ++it) {
      if (present[it - items.begin()] && p == *it) return true;
    }
    return false;
  }

  parlay::sequence<bool> contains(const parlay::sequence<objT> &points) const {
    parlay::sequence<bool> ret(points.size(), false);
    if (empty()) return ret;
    for (int i = 0; i < (int)points.size(); i++) {
      ret[i] = contains(points[i]);
    }
    return ret;
  }

  FoundPoint find(objT &p) const {
    if (!empty()) {
      auto node = nodes;
      nodeT *parent = nullptr;
      nodeT *gparent = nullptr;
      while (!node->isLeaf()) {
        gparent = parent;
        parent = node;

        // move down
        auto split_dim = node->getSplitDimension();
        auto split_val = node->getSplitValue();
        if (p.at(split_dim) < split_val)
          node = node->getLeft();
        else
          node = node->getRight();
      }

      // DEBUG_MSG(
      //"Find Checking: " << (FoundPoint){node - nodes, parent - nodes, gparent - nodes, -99});
      // [node] is a leaf -> check if it contains [p]
      for (auto it = node->getStartValue(); it != node->getEndValue(); ++it) {
        auto subtree_idx = (it - node->getStartValue());
        assert(subtree_idx >= 0);
        assert((size_t)subtree_idx < node->getValues().size());
        if (present[it - items.begin()] && p == *it) {
          return {node - nodes, parent - nodes, gparent - nodes, (int)subtree_idx};
        }
      }
    }
    return {FoundPoint::NOT_FOUND, FoundPoint::NOT_FOUND, FoundPoint::NOT_FOUND, -1};
  }

  parlay::sequence<objT> orthogonalQuery(objT &qMin, objT &qMax) {
    parlay::sequence<objT> ret;
    if (!empty()) {
      nodes[0].orthogonalQuery(qMin, qMax, items.begin(), present, ret);
    }
    return ret;
  }

#if (DUAL_KNN_MODE != DKNN_ARRAY)
  void updateDualDist(const parlay::slice<knnBuf::buffer<const pointT *> *,
                                          knnBuf::buffer<const pointT *> *> &buf_slice) {
    nodes[0].updateDualDist(buf_slice, items.begin());
  }
#else
  void updateDualDist(const parlay::slice<knnBuf::buffer<const pointT *> *,
                                          knnBuf::buffer<const pointT *> *> &buf_slice,
                      parlay::sequence<double> &dualKnnDists) {
    nodes[0].updateDualDist(buf_slice, items.begin(), nodes, dualKnnDists);
  }
#endif

  template <bool update, bool recurse_sibling>
  void knnSinglePoint(objT &p, knnBuf::buffer<pointT *> &buf) {
    nodes[0].template knnHelper<update, recurse_sibling>(
        pointT(p.coords()), items.begin(), present, buf);
    buf.keepK();  // TODO: could cause problems in logtree when running on nearly depleted
                  // subtree
  }

  template <bool set_res, bool update, bool recurse_sibling>
  void knnSinglePoint(
      objT &p,
      int i,
      parlay::slice<knnBuf::elem<pointT *> *, knnBuf::elem<pointT *> *> &out,
      parlay::slice<pointT **, pointT **> &res,
      int k,
      bool preload) {
    auto buf = knnBuf::buffer<pointT *>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
    if (preload) buf.ptr = k;
    knnSinglePoint<update, recurse_sibling>(p, buf);

    if (set_res) {
      for (int j = 0; j < k; j++) {
        res[i * k + j] = buf[j].entry;
      }
    }
  }

  template <bool update, bool recurse_sibling>
  void knn(const parlay::sequence<objT> &queries,
           parlay::slice<knnBuf::buffer<const pointT *> *, knnBuf::buffer<const pointT *> *> &bufs)
      const {
    assert(bufs.size() == queries.size());

    if (parallel) {
      parlay::parallel_for(0, queries.size(), [&](size_t i) {
        knnSinglePoint<update, recurse_sibling>(queries[i], bufs[i]);
      });
    } else {
      for (size_t i = 0; i < queries.size(); i++)
        knnSinglePoint<update, recurse_sibling>(queries[i], bufs[i]);
    }
  }

  template <bool set_res, bool update, bool recurse_sibling>
  void knn(parlay::sequence<objT> &queries,
           parlay::slice<knnBuf::elem<pointT *> *, knnBuf::elem<pointT *> *> &out,
           parlay::slice<pointT **, pointT **> &res,
           int k,
           bool preload = false) {
    assert(res.size() == k * queries.size());
    assert(out.size() == 2 * k * queries.size());

    if (parallel) {
      parlay::parallel_for(0, queries.size(), [&](size_t i) {
        knnSinglePoint<set_res, update, recurse_sibling>(queries[i], i, out, res, k, preload);
      });
    } else {
      for (size_t i = 0; i < queries.size(); i++)
        knnSinglePoint<set_res, update, recurse_sibling>(queries[i], i, out, res, k, preload);
    }
  }

  template <bool update, bool recurse_sibling>
  parlay::sequence<const pointT *> knn2(__attribute__((unused))
                                        const parlay::sequence<objT> &queries,
                                        __attribute__((unused)) int k) const {
    throw std::runtime_error("Called knn2 on the wrong tree type!");
  }
  template <bool update, bool recurse_sibling>
  parlay::sequence<const pointT *> knn3(__attribute__((unused))
                                        const parlay::sequence<objT> &queries,
                                        __attribute__((unused)) int k) const {
    throw std::runtime_error("Called knn3 on the wrong tree type!");
  }

  template <bool update = false, bool recurse_sibling = false>
  parlay::sequence<pointT *> knn(parlay::sequence<objT> &queries, int k) {
    parlay::sequence<pointT *> res(k * queries.size());
    parlay::sequence<knnBuf::elem<pointT *>> out(2 * k * queries.size());

    auto res_slice = res.head(res.size());
    auto out_slice = out.head(out.size());
    knn<true, update, recurse_sibling>(queries, out_slice, res_slice, k);
    return res;
  }

  // Dual knn stuff
  template <int _dim, class _objT, bool _parallel, bool _coarsen>
  friend parlay::sequence<const pargeo::point<_dim> *> dualKnn(
      parlay::sequence<_objT> &queries,
      const KdTree<_dim, _objT, _parallel, _coarsen> &rTree,
      int k);

  template <int _NUM_TREES,
            int _BUFFER_LOG2_SIZE,
            int _dim,
            class _objT,
            bool _parallel,
            bool _coarsen>
  friend parlay::sequence<const pargeo::point<_dim> *> dualKnn(
      parlay::sequence<_objT> &queries,
      const LogTree<_NUM_TREES, _BUFFER_LOG2_SIZE, _dim, _objT, _parallel, _coarsen> &rTree,
      int k);

#if (DUAL_KNN_MODE == DKNN_ARRAY)
  template <int _dim, class _objT, bool _parallel, bool _coarsenq, bool _coarsenr>
  friend void DualKnnHelper(kdNode<_dim, _objT, _parallel> *Q,
                            const kdNode<_dim, _objT, _parallel> *R,
                            const KdTree<_dim, _objT, _parallel, _coarsenq> &qTree,
                            parlay::sequence<double> &dualKnnDists,
                            const KdTree<_dim, _objT, _parallel, _coarsenr> &rTree,
                            parlay::slice<knnBuf::buffer<const point<_dim> *> *,
                                          knnBuf::buffer<const point<_dim> *> *> &bufs);
#else
  template <int _dim, class _objT, bool _parallel, bool _coarsenq, bool _coarsenr>
  friend void DualKnnHelper(kdNode<_dim, _objT, _parallel> *Q,
                            const kdNode<_dim, _objT, _parallel> *R,
                            const KdTree<_dim, _objT, _parallel, _coarsenq> &qTree,
                            const KdTree<_dim, _objT, _parallel, _coarsenr> &rTree,
                            parlay::slice<knnBuf::buffer<const pargeo::point<_dim> *> *,
                                          knnBuf::buffer<const pargeo::point<_dim> *> *> &bufs);
#endif

  parlay::sequence<const pointT *> dualKnnBase(const KdTree &queryTree, int k) const {
    parlay::sequence<const pointT *> res(k * queryTree.size());
    parlay::sequence<knnBuf::elem<const pointT *>> out(2 * k * queryTree.size());

    // allocate knn buffers
    parlay::sequence<knnBuf::buffer<const pointT *>> bufs(queryTree.size());
    if (parallel) {
      parlay::parallel_for(0, bufs.size(), [&](size_t i) {
        bufs[i] = knnBuf::buffer<const pointT *>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
      });
    } else {
      for (size_t i = 0; i < bufs.size(); i++) {
        bufs[i] = knnBuf::buffer<const pointT *>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
      }
    }

    // call dual knn helper
    auto buf_slice = bufs.cut(0, bufs.size());
#if (DUAL_KNN_MODE == DKNN_ARRAY)
    parlay::sequence<double> dualKnnDists(queryTree.num_nodes(),
                                          std::numeric_limits<double>::max());
    DualKnnHelper(queryTree.unsafe_root(), root(), queryTree, dualKnnDists, *this, buf_slice);
#else
    DualKnnHelper(queryTree.unsafe_root(), root(), queryTree, *this, buf_slice);
#endif

    // build result
    if (parallel) {
      parlay::parallel_for(0, queryTree.size(), [&](size_t i) {
        bufs[i].keepK();
        for (int j = 0; j < k; j++) {
          res[i * k + j] = bufs[i][j].entry;
        }
      });
    } else {
      for (size_t i = 0; i < queryTree.size(); i++) {
        bufs[i].keepK();
        for (int j = 0; j < k; j++) {
          res[i * k + j] = bufs[i][j].entry;
        }
      }
    }

    return res;
  }

  bool empty() const { return cur_size == 0; }
  size_t get_build_size() const { return build_size; }
  size_t size() const { return cur_size; }
  size_t capacity() const { return max_size; }
  size_t num_nodes() const { return 2 * max_size - 1; }
  auto node_idx(const nodeT *n) const {
    assert(n >= nodes);
    assert(n < nodes + num_nodes());
    return n - nodes;
  }

  // DELETE ========================================================================================
  // TODO: deduplicate points in bulk erase
  // Erase-By-Find: Erase the provided [FoundPoint] object. ----------------------------------------

  // TODO: this probably needs to propagate bbox recalculations all the way up the tree
  template <bool log_tree>
  void erase(const FoundPoint &found_point) {
    assert(found_point.idx != FoundPoint::NOT_FOUND);

    auto node = nodes + found_point.idx;
    auto parent = nodes + found_point.parent_idx;
    auto gparent = nodes + found_point.gparent_idx;

    assert(node->isLeaf());
    assert((parent == nodes) || (gparent != nullptr));  // either child of root or has grandparent

    // mark point as deleted
    present[found_point.point_idx + (node->getStartValue() - items.begin())] = false;
    node->removePoints(1);
    cur_size -= 1;

    // remove node if needed
    if (node->countPoints() == 0) {
      if (gparent != nullptr) {
        auto node_sibling = (parent->getLeft() == node) ? parent->getRight() : parent->getLeft();
        /* cut the parent out
         *         gp                   gp
         *        /  \                 /  \
         *       T    p      ==>      T    \
         *           / \                    \
         *          ns  n                   ns
         */
        if (gparent->getLeft() == parent) {
          gparent->setLeft(node_sibling);
        } else {
          gparent->setRight(node_sibling);
        }
        assert(!gparent->isLeaf());
        gparent->recomputeBoundingBox();
        // if (node_sibling) parents[node_sibling - nodes] = gparent;
      } else {
        if (!log_tree) {
          // only need to update structure if we're still gonna use it. -> in the case of a log
          // tree, if the entire side subtree of the root has been deleted, it will definitely be
          // moved up/down.
          if (parent->getLeft() == node)
            parent->setLeft(nullptr);
          else
            parent->setRight(nullptr);
          assert(!parent->isLeaf());
          parent->recomputeBoundingBox();
        }
      }
    } else {  // have to recompute bbox of node after erasing a point
      assert(node->isLeaf());
      node->recomputeBoundingBoxLeaf(items.begin(), present);
    }
  }

  // Basic Erase-By-Point: Directly erase the passed-in points in a single phase. ------------------
  template <bool log_tree = false>
  bool erase(objT &point) {
    const auto &found_node = find(point);
    if (found_node.idx != FoundPoint::NOT_FOUND) {
      // Found!
      erase<log_tree>(found_node);
      return true;
    }
    return false;
  }

  /*!
   * Delete input points (if they exist) from the tree.
   * @param points the list of points to be deleted.
   */
  template <bool log_tree>
  void erase(parlay::slice<objT *, objT *> points) {
    for (auto &p : points)
      erase<log_tree>(p);
  }

  template <bool log_tree>
  void erase(parlay::sequence<objT> &points) {
    for (auto &p : points)
      erase<log_tree>(p);
  }

  // Serial Erase Functions ------------------------------------------------------------------------
  //  - Provide an unoptimized, serial erasure operation
  //[[deprecated("basically the same as serial erase, above")]] parlay::sequence<FoundPoint>
  // serial_erase_prepare(const parlay::sequence<objT> &points) const {
  // parlay::sequence<FoundPoint> ret(
  // points.size(), {FoundPoint::NOT_FOUND, FoundPoint::NOT_FOUND, FoundPoint::NOT_FOUND});
  // if (empty()) return ret;
  // for (int i = 0; i < (int)points.size(); i++) {
  // ret[i] = find(points[i]);
  //}
  // return ret;
  //}

  /*!
   * Delete the input nodes from the tree. The input should be directly passed from
   * [serial_erase_prepare].
   * TODO: Currently assumes that input is valid - need to refactor FoundPoint to be opaque.
   * @param found_nodes the list of nodes to be deleted.
   */
  // template <bool destroy = false>
  //[[deprecated("serial deletion is a single phase; this doesn't work anymore")]] void
  // serial_erase_commit(parlay::slice<FoundPoint *, FoundPoint *> found_nodes) {
  // for (const auto &found_node : found_nodes) {
  //// have to recompute parents because the erasures are being done serially
  // auto idx = found_node.idx;
  // auto parent_idx = found_node.parent_idx;    // parents[idx] - nodes;
  // auto gparent_idx = found_node.gparent_idx;  // parents[parent_idx] - nodes;
  // erase<destroy>({idx, parent_idx, gparent_idx});
  //}
  //}

  // Bulk Erase Functions --------------------------------------------------------------------------
  /*!
   * Delete the input nodes from the tree. The input should be directly passed from
   * [bulk_erase].
   */
  // TODO: think about whether you want objT or objT* in points.
  nodeT *bulk_erase_leaf(nodeT *node, parlay::slice<objT *, objT *> points, size_t &num_removed) {
    assert(node->isLeaf());
    // remove all the points we can
    num_removed = 0;
#ifdef ERASE_SEARCH_TIMES
    timer t;
#endif
    for (const auto &pt_to_del : points) {
      for (auto it = node->getStartValue(); it != node->getEndValue(); ++it) {
        if (present[it - items.begin()] && pt_to_del == *it) {
          present[it - items.begin()] = false;
          num_removed++;
          break;
        }
      }
    }
#ifdef ERASE_SEARCH_TIMES
    total_search_time += t.get_next();
#endif

    node->removePoints(num_removed);  // mark the points as removed
#ifdef ERASE_SEARCH_TIMES
    total_bbox_time += t.get_next();
#endif

    // check if we can delete the leaf
    if (node->countPoints() == 0) {
      return nullptr;  // deleted all points -> delete the leaf
    } else {
      if (num_removed > 0) {
        node->recomputeBoundingBoxLeaf(items.begin(), present);
      }
      return node;  // didn't delete all the points -> don't delete this leaf
    }
  }

  nodeT *bulk_erase_helper(nodeT *node, parlay::slice<objT *, objT *> points, size_t &num_removed) {
    if (node->isLeaf()) {
#ifdef ERASE_SEARCH_TIMES
      timer t;
#endif
      auto r = bulk_erase_leaf(node, points, num_removed);
#ifdef ERASE_SEARCH_TIMES
      total_leaf_time += t.get_next();
#endif
      return r;
    } else {
      auto right_start = serialPartition(points, node->getSplitDimension(), node->getSplitValue());

      // recurse on the two halves
      size_t num_removed_left, num_removed_right;
      auto new_left =
          bulk_erase_helper(node->getLeft(), points.cut(0, right_start), num_removed_left);
      auto new_right = bulk_erase_helper(
          node->getRight(), points.cut(right_start, points.size()), num_removed_right);

      num_removed = num_removed_left + num_removed_right;
      // need to deal with root
      if (new_left != nullptr && new_right != nullptr) {
        // neither child deleted -> i'm not deleted
        // TODO: it's possible to further case based on whether left/right changed, but doesn't seem
        // worth
        node->setLeft(new_left);
        node->setRight(new_right);
        assert(!node->isLeaf());
        node->recomputeBoundingBox();
        // parents[new_left - nodes] = node;
        // parents[new_right - nodes] = node;
        return node;
      } else if (new_left == nullptr && new_right == nullptr) {
        // both children deleted -> delete me
        return nullptr;
      } else {
        // one child deleted -> delete me, but not my whole subtree
        return (new_left == nullptr) ? new_right : new_left;
      }
    }
  }

  nodeT *bulk_erase_helper_parallel(nodeT *node,
                                    parlay::slice<objT *, objT *> points,
                                    parlay::slice<bool *, bool *> flags,
                                    size_t &num_removed) {
    assert(parallel);
    if (!eraseInParallel(points.size())) {  // just fallback to serial
      return bulk_erase_helper(node, points, num_removed);
    }

#ifdef PRINT_KDTREE_TIMINGS
    bool print_timings = (points.size() == 2000000);
    auto mtime = [&](const std::string &s) {
      if (print_timings) this->mark_time(s);
    };
#endif

    if (node->isLeaf()) {
      return bulk_erase_leaf(node, points, num_removed);
    } else {
#ifdef PRINT_KDTREE_TIMINGS
      mtime("Start partition: " + std::to_string(points.size()));
#endif
      auto right_start =
          parallelPartition(points, flags, node->getSplitDimension(), node->getSplitValue());
#ifdef PRINT_KDTREE_TIMINGS
      mtime("Finish partition: " + std::to_string(points.size()));
#endif

      // recurse on the two halves
      nodeT *new_left, *new_right;
      size_t num_removed_left, num_removed_right;
      parlay::par_do(
          [&]() {
            new_left = bulk_erase_helper_parallel(node->getLeft(),
                                                  points.cut(0, right_start),
                                                  flags.cut(0, right_start),
                                                  num_removed_left);
          },
          [&]() {
            new_right = bulk_erase_helper_parallel(node->getRight(),
                                                   points.cut(right_start, points.size()),
                                                   flags.cut(right_start, points.size()),
                                                   num_removed_right);
          });
      num_removed = num_removed_left + num_removed_right;
#ifdef PRINT_KDTREE_TIMINGS
      mtime("Finish Recursion: " + std::to_string(points.size()));
#endif

      // need to deal with root
      if (new_left != nullptr && new_right != nullptr) {
        // neither child deleted -> i'm not deleted
        // TODO: it's possible to further case based on whether left/right changed, but doesn't seem
        // worth
        // TODO: can parallelize with atomics if we want
        node->setLeft(new_left);
        node->setRight(new_right);
        assert(!node->isLeaf());
        node->recomputeBoundingBox();
        // parents[new_left - nodes] = node;
        // parents[new_right - nodes] = node;
        return node;
      } else if (new_left == nullptr && new_right == nullptr) {
        // both children deleted -> delete me
        return nullptr;
      } else {
        // one child deleted -> delete me, but not my whole subtree
        return (new_left == nullptr) ? new_right : new_left;
      }
    }
  }

  // Have to make a copy of the points array so that we can move them around.
#ifdef ALL_USE_BLOOM
  void bulk_erase(const parlay::sequence<objT> &points_in)
#else
  void bulk_erase(parlay::sequence<objT> &points)
#endif
  {
    if (empty()) return;
      // One scenario where this could create an invalid tree: if the entire left/right subtree
      // of the root is being deleted, this function will not update root. BUT we don't care about
      // this scenario, because this necessarily implies that the tree is being rebuilt.
#ifdef ALL_USE_BLOOM
    auto points = bloom_filter.filter(points_in);
#endif
    size_t num_removed;
    if (parallel) {
#ifdef PRINT_KDTREE_TIMINGS
      this->mark_time("Erase Start");
#endif
      auto flags = parlay::sequence<bool>(points.size());
      bulk_erase_helper_parallel(
          nodes, points.cut(0, points.size()), flags.cut(0, points.size()), num_removed);
#ifdef PRINT_KDTREE_TIMINGS
      this->mark_time("Erase Finish");
#endif
    } else {
      bulk_erase_helper(nodes, points.cut(0, points.size()), num_removed);
    }
    cur_size -= num_removed;
  }

  // DEBUG --------------------------------------------
  /*!
   * Verify that the tree is balanced - the right child has either the same number of points as the
   * left, or exactly one more.
   */
  bool verify() const {
    nodes[0].verify();
    return true;
  }

  std::pair<parlay::slice<const objT *, const objT *>, parlay::slice<const bool *, const bool *>>
  getItems() const {
    return {items.cut(0, build_size), present.cut(0, build_size)};
  }

  /*!
   * Get direct pointer to the root node of the tree.
   * @return the root node (pointer to node array).
   */
  const nodeT *root() const { return nodes; }
  nodeT *unsafe_root() const { return nodes; }

  // https://stackoverflow.com/questions/36802354/print-binary-tree-in-a-pretty-way-using-c
  void print(const std::string &prefix, const nodeT *node, bool isLeft) const {
    if (node != nullptr) {
      std::cout << prefix;
      std::cout << (isLeft ? "├──" : "└──");
      node->print(getNodeValueIdx(node));

      // print next level
      print(prefix + (isLeft ? "│   " : "    "), node->getLeft(), true);
      print(prefix + (isLeft ? "│   " : "    "), node->getRight(), false);
    }
  }
  void print() const { print(std::string(""), nodes, false); }
};

#endif  // KDTREE_H
