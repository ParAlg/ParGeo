#ifndef KDNODE_H
#define KDNODE_H

#include <atomic>
#include "parlay/parallel.h"
#include "parlay/sequence.h"
//#include "common/geometry.h"
#include "pargeo/point.h"

#include "macro.h"
#include "knnbuffer.h"
#include "box.h"

template <int dim, class objT, bool parallel, bool coarsen>
class KdTree;

// make dim intrinsic to objt todo
template <int dim, class objT, bool parallel>
class kdNode {
  typedef int intT;
  typedef double floatT;
  typedef pargeo::point<dim> pointT;
  typedef kdNode<dim, objT, parallel> nodeT;

  // TODO: split leaf/non-leaf node data
  // leaf nodes
  int num_points;
  // non-leaf node
  int split_dimension;  // TODO: think about making this dynamic, instead of alternating
  floatT split_value;

  // all nodes
  parlay::slice<objT *, objT *> subtree_items;  // TODO: make these const pointers
  pointT pMin, pMax;

  // Node pointers
  nodeT *left;
  nodeT *right;

  // dual knn distances (only used for queries)
#if (DUAL_KNN_MODE == DKNN_ATOMIC_LEAF)
  std::atomic<double> dualKnnDist;
#elif (DUAL_KNN_MODE == DKNN_NONATOMIC_LEAF)
  double dualKnnDist;
#endif

  bool computeRangeQueryInParallel() const {
    return false;
    if (!left || !right) return false;  // only one child
    return subtree_items.size() >= RANGEQUERY_BASE_CASE;
  }
  bool computeBoundingBoxInParallel() const {
    if (!left || !right) return false;  // only one child
    return subtree_items.size() >= BOUNDINGBOX_BASE_CASE;
  }
  bool dualKnnRecurseInParallel() const {
    if (!left || !right) return false;  // only one child
    return subtree_items.size() >= DUALKNN_BASE_CASE;
  }

 public:
  // non-leaf
  kdNode(int split_dimension_, floatT split_value_, parlay::slice<objT *, objT *> subtree_items_)
      : num_points(-1),
        split_dimension(split_dimension_),
        split_value(split_value_),
        subtree_items(std::move(subtree_items_)),
        left(nullptr),
        right(nullptr)
#if (DUAL_KNN_MODE != DKNN_ARRAY)
        ,
        dualKnnDist(std::numeric_limits<double>::max())
#endif
  {
  }
  // leaf
  kdNode(parlay::slice<objT *, objT *> subtree_items_) : kdNode(-1, 0, subtree_items_) {
    assert(subtree_items.size() > 0);
    num_points = subtree_items.size();
    pMin = pointT(subtree_items[0].coords());
    pMax = pointT(subtree_items[0].coords());
    for (auto &pt : subtree_items) {
      pMin.minCoords(pt.coords());
      pMax.maxCoords(pt.coords());
    }
  }

  // Modifiers
#if (DUAL_KNN_MODE != DKNN_ARRAY)
  void updateDualDist(const parlay::slice<knnBuf::buffer<const pointT *> *,
                                          knnBuf::buffer<const pointT *> *> &buf_slice,
                      const objT *tree_start) {
#else
  void updateDualDist(const parlay::slice<knnBuf::buffer<const pointT *> *,
                                          knnBuf::buffer<const pointT *> *> &buf_slice,
                      const objT *tree_start,
                      const kdNode<dim, objT, parallel> *nodes,
                      parlay::sequence<double> &dualKnnDists) {
#endif
    double new_rad = 0;
    if (isLeaf() && buf_slice[subtree_items.begin() - tree_start].hasK()) {
      for (auto it = subtree_items.begin(); it != subtree_items.end(); ++it) {
        auto idx = it - tree_start;
        new_rad = std::max(new_rad, buf_slice[idx].keepK().cost);
      }
    } else {
      auto update_node = [&](kdNode<dim, objT, parallel> *n) {
#if (DUAL_KNN_MODE != DKNN_ARRAY)
        n->updateDualDist(buf_slice, tree_start);
#else
        n->updateDualDist(buf_slice, tree_start, nodes, dualKnnDists);
#endif
      };

      auto node_dist = [&](const kdNode<dim, objT, parallel> *n) {
#if (DUAL_KNN_MODE == DKNN_ATOMIC_LEAF)
        return n->dualKnnDist.load();
#elif (DUAL_KNN_MODE == DKNN_NONATOMIC_LEAF)
        return n->dualKnnDist;
#else
      return dualKnnDists[n - nodes];
#endif
      };

      if (parallel && computeBoundingBoxInParallel()) {
        parlay::par_do([&]() { update_node(left); }, [&]() { update_node(right); });
        new_rad = std::max(node_dist(left), node_dist(right));
      } else {
        if (left) {
          update_node(left);
          new_rad = std::max(new_rad, node_dist(left));
        }
        if (right) {
          update_node(right);
          new_rad = std::max(new_rad, node_dist(right));
        }
      }
    }
#if (DUAL_KNN_MODE != DKNN_ARRAY)
    update_dual_knn_dist(new_rad);
#else
    dualKnnDists[this - nodes] = new_rad;
#endif
  }

#if (DUAL_KNN_MODE != DKNN_ARRAY)
  inline void update_dual_knn_dist(const double new_val) {
#if (DUAL_KNN_MODE == DKNN_ATOMIC_LEAF)
    // https://stackoverflow.com/questions/44800510/updating-maximum-value-atomically
    double cur_val = dualKnnDist;
    while (cur_val > new_val) {
      dualKnnDist.compare_exchange_weak(cur_val, new_val, std::memory_order_relaxed);
    }
#else
    dualKnnDist = new_val;
#endif
  }
#endif

  void removePoints(int num) {
    assert(isLeaf());
    assert(num >= 0);
    assert(num <= num_points);
    num_points -= num;
  }
  void setLeft(nodeT *p) { left = p; }
  void setRight(nodeT *p) { right = p; }
  void setEmpty() { split_dimension = -2; }

  void recomputeBoundingBox() {
    // assumes child bounding boxes are computed
    if (left && right) {
      pMin = pointT(left->getMin().coords());
      pMax = pointT(left->getMax().coords());
      pMin.minCoords(right->getMin().coords());
      pMax.maxCoords(right->getMax().coords());
    } else if (left) {
      pMin = pointT(left->getMin().coords());
      pMax = pointT(left->getMax().coords());
    } else if (right) {
      pMin = pointT(right->getMin().coords());
      pMax = pointT(right->getMax().coords());
    } else {
      assert(isLeaf());
    }
  }

  void recomputeBoundingBoxLeaf(const objT *tree_start, const parlay::sequence<bool> &present) {
    assert(isLeaf());
    bool first = true;
    for (auto it = subtree_items.begin(); it != subtree_items.end(); ++it) {
      if (present[it - tree_start]) {
        if (first) {
          pMin = pointT(it->coords());
          pMax = pointT(it->coords());
          first = false;
        } else {
          pMin.minCoords(it->coords());
          pMax.maxCoords(it->coords());
        }
      }
    }
  }

  // recompute for entire subtree
  void recomputeBoundingBoxSubtree() {
    if (parallel && computeBoundingBoxInParallel()) {
      parlay::par_do([&]() { left->recomputeBoundingBoxSubtree(); },
                     [&]() { right->recomputeBoundingBoxSubtree(); });
    } else {
      if (left) left->recomputeBoundingBoxSubtree();
      if (right) right->recomputeBoundingBoxSubtree();
    }
    recomputeBoundingBox();
  }

  // Getters
  pointT &getMin() { return pMin; }
  pointT &getMax() { return pMax; }
  nodeT *getLeft() { return left; }
  nodeT *getRight() { return right; }
  bool isLeaf() const { return (left == nullptr) && (right == nullptr); }
  bool isEmpty() const { return split_dimension == -2; }

  int countPoints() const {
    if (isLeaf()) return num_points;
    return left->countPoints() + right->countPoints();
  }

  const auto &getValues() const {
    assert(isLeaf());
    return subtree_items;
  }

  const objT *getStartValue() const { return subtree_items.begin(); }
  const objT *getEndValue() const { return subtree_items.end(); }

  int getSplitDimension() const {
    assert(!isLeaf());
    return split_dimension;
  }
  floatT getSplitValue() const {
    assert(!isLeaf());
    return split_value;
  }

  // Query
  // MOVED - [contains] is performed in [KdTree]
  // bool contains(const objT &p, const parlay::sequence<bool> &present) const {
  // if (isLeaf()) {
  // for (const auto &pt : subtree_items) {
  // if (pt == p) return true;
  //}
  // return false;
  //}

  //// internal (non-leaf) node
  // if (p.coordinate(split_dimension) < split_value) return left->contains(p);
  // return right->contains(p);
  //}

  // TODO: can probably make this recurse more intelligently if we precompute return sizes
  void orthogonalQuery(objT &qMin,
                       objT &qMax,
                       objT *tree_start,
                       parlay::sequence<bool> &present,
                       parlay::sequence<objT> &ret) const {
    auto cmp = boxCompare(qMin, qMax, pMin, pMax);
    if (cmp == BOX_EXCLUDE) {
      return;
    } else if (cmp == BOX_INCLUDE) {  // query box contains node box -> take all the points
      // allocate space for new points
      auto orig_ret_size = ret.size();
      ret.resize(ret.size() + subtree_items.size());  // TODO: precompute the exact size

      // compute [present] subarray
      auto start = getStartValue() - tree_start;
      auto end = getEndValue() - tree_start;
      assert(end > start);
      assert(subtree_items.size() == (size_t)(end - start));

      auto num_added = parlay::pack_into(
          subtree_items, present.cut(start, end), ret.cut(orig_ret_size, ret.size()));

      // resize
      ret.resize(orig_ret_size + num_added);
    } else {
      assert(cmp == BOX_OVERLAP);
      if (isLeaf()) {
        // TODO: maybe do this more intelligently? (precompute and/or parallelize)
        for (auto it = subtree_items.begin(); it != subtree_items.end(); ++it) {
          if (present[it - tree_start] && itemInBox(qMin, qMax, it)) {
            ret.push_back(it);
          }
        }
      } else if (parallel && computeRangeQueryInParallel()) {
        assert(left);
        assert(right);

        parlay::sequence<objT> right_ret;

        parlay::par_do(
            [&]() { left->orthogonalQuery(qMin, qMax, tree_start, present, ret); },
            [&]() { right->orthogonalQuery(qMin, qMax, tree_start, present, right_ret); });

        // put right_ret into ret
        ret.insert(ret.begin() + ret.size(), right_ret.begin(), right_ret.end());
      } else {
        if (left) left->orthogonalQuery(qMin, qMax, tree_start, present, ret);
        if (right) right->orthogonalQuery(qMin, qMax, tree_start, present, ret);
      }
    }
  }

  void knnAddToBuffer(const pointT &q,
                      const objT *tree_start,
                      const parlay::sequence<bool> &present,
                      knnBuf::buffer<const pointT *> &out,
                      double radius = std::numeric_limits<double>::max()) const {
    // TODO: maybe parallelize?
    auto start = getStartValue() - tree_start;
    [[maybe_unused]] auto end = getEndValue() - tree_start;
    assert(end > start);
    assert(subtree_items.size() == (size_t)(end - start));

    for (size_t i = 0; i < subtree_items.size(); i++) {
      if (present[start + i]) {  // point isn't deleted
        auto dist = q.dist(subtree_items[i]);
        if (dist <= radius) {  // point within radius of interest
          const pointT *item_ptr = subtree_items.begin() + i;
          out.insert(knnBuf::elem(dist, item_ptr));
        }
      }
    }
  }

  template <bool update>
  void knnPrune(const pointT &q,
                const objT *tree_start,
                const parlay::sequence<bool> &present,
                double &radius,
                pointT &qMin,
                pointT &qMax,
                knnBuf::buffer<const pointT *> &out) const {
    if (update) {
      // compute current radius
      auto tmp = out.keepK();
      auto new_radius = tmp.cost;

      // update the query box if necessary
      if (new_radius < radius) {
        radius = new_radius;
        // create box based on radius
        for (int i = 0; i < dim; i++) {
          qMin[i] = q.at(i) - radius;
          qMax[i] = q.at(i) + radius;
        }
      }
    }

    // search only the intersection of the subtree with the radius-box
    auto cmp = boxCompare(qMin, qMax, pMin, pMax);
    switch (cmp) {
      case BOX_EXCLUDE: {
        return;
      }
      case BOX_INCLUDE: {
        knnAddToBuffer(q, tree_start, present, out, radius);
        break;
      }
      case BOX_OVERLAP: {
        if (isLeaf()) {
          knnAddToBuffer(q, tree_start, present, out, radius);
        } else {
          left->knnPrune<update>(q, tree_start, present, radius, qMin, qMax, out);
          right->knnPrune<update>(q, tree_start, present, radius, qMin, qMax, out);
        }
        break;
      }
    }
  }

  // Taken with modifications from:
  // https://github.mit.edu/yiqiuw/pargeo/blob/master/knnSearch/kdTree/kdtKnn.h#L365
  template <bool update, bool recurse_sibling>
  void knnHelper(const pointT &q,
                 const objT *tree_start,
                 const parlay::sequence<bool> &present,
                 knnBuf::buffer<const pointT *> &out) const {
    // first, find the leaf
    nodeT *other_child;
    if (isLeaf()) {
      knnAddToBuffer(q, tree_start, present, out);
      return;  // base case
    } else {
      if (q.at(split_dimension) < split_value) {
        // TODO: hint to compiler that [left] will pretty much never be null
        if (left) left->knnHelper<update, recurse_sibling>(q, tree_start, present, out);
        other_child = right;
      } else {
        if (right) right->knnHelper<update, recurse_sibling>(q, tree_start, present, out);
        other_child = left;
      }
    }

    // now, check alternate children with aggressive pruning
    if (!out.hasK()) {
      // try finding knn on other child
      if (recurse_sibling) {
        other_child->knnHelper<update, recurse_sibling>(q, tree_start, present, out);
      } else {
        other_child->knnAddToBuffer(q, tree_start, present, out);
      }
    } else {
      double radius = std::numeric_limits<double>::max();
      pointT qMin, qMax;

      if (!update) {
        // compute current radius
        auto tmp = out.keepK();
        auto new_radius = tmp.cost;

        // update the query box if necessary
        if (new_radius < radius) {
          radius = new_radius;
          // create box based on radius
          for (int i = 0; i < dim; i++) {
            qMin[i] = q.at(i) - radius;
            qMax[i] = q.at(i) + radius;
          }
        } else {
          assert(false);
        }
      }

      other_child->knnPrune<update>(q, tree_start, present, radius, qMin, qMax, out);
    }
  }

  // Dual knn stuff
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

  // Debug
  // TODO: make this also check that points are on the right side of a split
  int verify() const {
    if (isLeaf()) return countPoints();

    auto left_points = left->verify();
    auto right_points = right->verify();
    // the children must be equal sized, or right has one more item
    if (!((left_points == right_points) || (left_points + 1 == right_points)))
      throw std::runtime_error("Invalid tree!: (left#, right#) = (" + std::to_string(left_points) +
                               ", " + std::to_string(right_points) + ")");
    return left_points + right_points;
  }

  void print(const std::pair<int, int> &idx) const {
    std::stringstream ss;
    ss << "idx = [" << idx.first << ", " << idx.second << ")";
    if (isLeaf()) {
      std::cout << "leaf: { " << ss.str() << "; ";
      for (const auto &pt : subtree_items) {
        std::cout << "(";
        for (int i = 0; i < dim; i++) {
          std::cout << pt->at(i);
          if (i < dim - 1) std::cout << ", ";
        }
        std::cout << "), ";
      }
      std::cout << " }" << std::endl;
    } else {
      auto split_dim = getSplitDimension();
      auto split_val = getSplitValue();
      std::cout << "split: { " << ss.str() << "; dim = " << split_dim << "; val = " << split_val
                << "}" << std::endl;
    }
  }
};

// Helper function to find bounding box distances between nodes
template <int dim, class objT, bool parallel>
inline double KdNodeBoundingBoxDistance(const kdNode<dim, objT, parallel> *n1,
                                        const kdNode<dim, objT, parallel> *n2) {
  return BoundingBoxDistance(n1->getMin(), n1->getMax(), n2->getMin(), n2->getMax());
}

// For now, this is not cache-oblivious. Instead, it uses a simple, d-dimension generalizable
// approach.
template <int dim, class objT>
class Auxiliary {
  parlay::sequence<objT *> items;
  parlay::sequence<parlay::sequence<objT *>> Ph, Pv;
  parlay::sequence<objT *> sortedX, sortedY;

 public:
  Auxiliary(parlay::sequence<objT *> items_)
      : items(std::move(items_)) {}  // create internal copy of items
};

#endif  // KDNODE_H
