#ifndef DUAL_H
#define DUAL_H

#include "../cache-oblivious/cokdtree.h"
#include "../log-tree/logtree.h"

#ifdef PRINT_DKNN_TIMINGS
#include "common/get_time.h"
#endif

// Top-level wrappers for calling dual knn
template <int dim, class objT, bool parallel, bool coarsen>
parlay::sequence<const pargeo::point<dim> *> dualKnn(parlay::sequence<objT> &queries,
                                             const KdTree<dim, objT, parallel, coarsen> &rTree,
                                             int k) {
  // construct query tree
#ifdef PRINT_DKNN_TIMINGS
  timer t;
#endif
  CO_KdTree<dim, objT, parallel, coarsen> qTree((int)std::ceil(std::log2(queries.size())));
  qTree.build(std::move(queries));
#ifdef PRINT_DKNN_TIMINGS
  std::cout << "[DKNN] Query Tree Construction: " << t.get_next() << "\n";
#endif

  auto ret = rTree.dualKnnBase(qTree, k);  // call dual knn
  queries = std::move(qTree.items);        // move the query points back

  return ret;
}

template <int NUM_TREES,         // the number of static trees
          int BUFFER_LOG2_SIZE,  // the size of the (dynamic) buffer tree
          int dim,
          class objT,
          bool parallel,
          bool coarsen>
parlay::sequence<const pargeo::point<dim> *> dualKnn(
    parlay::sequence<objT> &queries,
    const LogTree<NUM_TREES, BUFFER_LOG2_SIZE, dim, objT, parallel, coarsen> &rTree,
    int k) {
  // construct query tree
#ifdef PRINT_DKNN_TIMINGS
  timer t;
#endif
  CO_KdTree<dim, objT, parallel, coarsen> qTree((int)std::ceil(std::log2(queries.size())));
  qTree.build(std::move(queries));
#ifdef PRINT_DKNN_TIMINGS
  std::cout << "[DKNN] Query Tree Construction: " << t.get_next() << "\n";
#endif

  auto ret = rTree.dualKnnBase(qTree, k);  // call dual knn
  queries = std::move(qTree.items);        // move the query points back

  return ret;
}

template <int dim, class objT, bool parallel, bool coarsenq, bool coarsenr>
#if (DUAL_KNN_MODE == DKNN_ARRAY)
void DualKnnHelper(kdNode<dim, objT, parallel> *Q,
                   const kdNode<dim, objT, parallel> *R,
                   const KdTree<dim, objT, parallel, coarsenq> &qTree,
                   parlay::sequence<double> &dualKnnDists,
                   const KdTree<dim, objT, parallel, coarsenr> &rTree,
                   parlay::slice<knnBuf::buffer<const point<dim> *> *,
                                 knnBuf::buffer<const point<dim> *> *> &bufs) {
#else
void DualKnnHelper(kdNode<dim, objT, parallel> *Q,
                   kdNode<dim, objT, parallel> *R,
                   KdTree<dim, objT, parallel, coarsenq> &qTree,
                   KdTree<dim, objT, parallel, coarsenr> &rTree,
                   parlay::slice<knnBuf::buffer<const pargeo::point<dim> *> *,
                                 knnBuf::buffer<const pargeo::point<dim> *> *> &bufs) {
#endif
  typedef kdNode<dim, objT, parallel> nodeT;
  // if (!Q || !R) return;
  assert(Q && R);

  auto recurse = [&](kdNode<dim, objT, parallel> *_Q, const kdNode<dim, objT, parallel> *_R) {
#if (DUAL_KNN_MODE == DKNN_ARRAY)
    DualKnnHelper(_Q, _R, qTree, dualKnnDists, rTree, bufs);
#else
    DualKnnHelper(_Q, _R, qTree, rTree, bufs);
#endif
  };

  // if either node has only a single child, just forward the call
  if (Q->getLeft() && !Q->getRight()) {
    recurse(Q->getLeft(), R);
    return;
  } else if (!Q->getLeft() && Q->getRight()) {
    recurse(Q->getRight(), R);
    return;
  } else if (R->getLeft() && !R->getRight()) {
    recurse(Q, R->getLeft());
    return;
  } else if (!R->getLeft() && R->getRight()) {
    recurse(Q, R->getRight());
    return;
  }

  // either a leaf or a node with 2 children
  assert(Q->isLeaf() || (Q->getLeft() && Q->getRight()));
  assert(R->isLeaf() || (R->getLeft() && R->getRight()));

  // order calls based on bbox distance (TODO: does this actually help?)
  // used below for the one-sided recursion cases
  auto one_sided_recurse =
      [&](bool recurseInParallel, nodeT *Q1, const nodeT *R1, nodeT *Q2, const nodeT *R2) {
        auto dist1 = KdNodeBoundingBoxDistance(Q1, R1);
        auto dist2 = KdNodeBoundingBoxDistance(Q2, R2);
        if (dist1 < dist2) {  // 1 before 2
          if (recurseInParallel) {
            parlay::par_do([&]() { recurse(Q1, R1); }, [&]() { recurse(Q2, R2); });
          } else {
            recurse(Q1, R1);
            recurse(Q2, R2);
          }
        } else {  // 2 before 1
          if (recurseInParallel) {
            parlay::par_do([&]() { recurse(Q2, R2); }, [&]() { recurse(Q1, R1); });
          } else {
            recurse(Q2, R2);
            recurse(Q1, R1);
          }
        }
      };

#if (DUAL_KNN_MODE == DKNN_ARRAY)
  if (KdNodeBoundingBoxDistance(Q, R) > dualKnnDists[qTree.node_idx(Q)]) {
#else
  if (KdNodeBoundingBoxDistance(Q, R) > Q->dualKnnDist) {
#endif
    // definitely no updates here
    return;
  } else if (Q->isLeaf() && R->isLeaf()) {  // reached leaves -> update!
    // get the indices of the items in Q
    const auto &q_items = qTree.items;
    auto q_start = Q->getStartValue() - q_items.begin();
    auto q_end = Q->getEndValue() - q_items.begin();
    assert(q_end > q_start);
    assert(q_start >= 0);
    assert(Q->subtree_items.size() == (size_t)(q_end - q_start));
    assert(qTree.size() == q_items.size());

    // leaves are small -> don't bother parallelizing
    double Q_new_radius = 0;
    for (auto q_idx = (size_t)q_start; q_idx < (size_t)q_end; q_idx++) {
      assert(q_idx >= 0);
      assert(q_idx < bufs.size());
      auto &q_out = bufs[q_idx];
      auto q_radius = q_out.hasK() ? q_out.keepK().cost : std::numeric_limits<double>::max();

      // relax in all the points in R
      R->knnAddToBuffer(q_items[q_idx], rTree.items.begin(), rTree.present, q_out, q_radius);

      // update the new radius
      auto q_new_radius = q_out.hasK() ? q_out.keepK().cost : std::numeric_limits<double>::max();
      Q_new_radius = std::max(Q_new_radius, q_new_radius);
    }

#if (DUAL_KNN_MODE == DKNN_ARRAY)
    dualKnnDists[qTree.node_idx(Q)] = Q_new_radius;
#else
    Q->update_dual_knn_dist(Q_new_radius);
#endif
  } else if (Q->isLeaf()) {
    // cannot recurse in parallel because Q is the same in both cases
    one_sided_recurse(false, Q, R->getLeft(), Q, R->getRight());
  } else if (R->isLeaf()) {
    one_sided_recurse(parallel && Q->dualKnnRecurseInParallel(), Q->getLeft(), R, Q->getRight(), R);

#if (DUAL_KNN_MODE == DKNN_ARRAY)
    dualKnnDists[qTree.node_idx(Q)] =
        std::max(dualKnnDists[qTree.node_idx(Q->left)], dualKnnDists[qTree.node_idx(Q->right)]);
#else
    Q->update_dual_knn_dist(std::max(Q->left->dualKnnDist, Q->right->dualKnnDist));
#endif
  } else {  // neither is leaf, all 4 recursive steps
    auto QlRl_dist = KdNodeBoundingBoxDistance(Q->getLeft(), R->getLeft());
    auto QlRr_dist = KdNodeBoundingBoxDistance(Q->getLeft(), R->getRight());
    // closer R child to Q->getLeft()
    auto Ql_R1 = (QlRl_dist < QlRr_dist) ? R->getLeft() : R->getRight();
    // further R child to Q->getLeft()
    auto Ql_R2 = (QlRl_dist < QlRr_dist) ? R->getRight() : R->getLeft();

    auto QrRl_dist = KdNodeBoundingBoxDistance(Q->getRight(), R->getLeft());
    auto QrRr_dist = KdNodeBoundingBoxDistance(Q->getRight(), R->getRight());
    auto Qr_R1 = (QrRl_dist < QrRr_dist) ? R->getLeft() : R->getRight();
    auto Qr_R2 = (QrRl_dist < QrRr_dist) ? R->getRight() : R->getLeft();

    if (parallel && (Q->dualKnnRecurseInParallel() || R->dualKnnRecurseInParallel())) {
      parlay::par_do([&]() { recurse(Q->getLeft(), Ql_R1); },
                     [&]() { recurse(Q->getRight(), Qr_R1); });
      parlay::par_do([&]() { recurse(Q->getLeft(), Ql_R2); },
                     [&]() { recurse(Q->getRight(), Qr_R2); });
    } else {
      recurse(Q->getLeft(), Ql_R1);
      recurse(Q->getRight(), Qr_R1);
      recurse(Q->getLeft(), Ql_R2);
      recurse(Q->getRight(), Qr_R2);
    }
#if (DUAL_KNN_MODE == DKNN_ARRAY)
    dualKnnDists[qTree.node_idx(Q)] =
        std::max(dualKnnDists[qTree.node_idx(Q->left)], dualKnnDists[qTree.node_idx(Q->right)]);
#else
    Q->update_dual_knn_dist(std::max(Q->left->dualKnnDist, Q->right->dualKnnDist));
#endif
  }
}

#endif  // DUAL_H
