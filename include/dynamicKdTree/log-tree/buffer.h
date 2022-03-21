#ifndef LOGTREE_BUFFER_H
#define LOGTREE_BUFFER_H

#include <parlay/parallel.h>
#include <parlay/sequence.h>

#include "../shared/macro.h"

template <int dim, class objT, bool parallel = false>
class alignas(64) LogTreeBuffer {
  typedef point<dim> pointT;
  // simple wrapper around parlay::sequence
  parlay::sequence<objT> items;
  parlay::sequence<bool> present;
  size_t cur_size;
  size_t insert_size;

 public:
  LogTreeBuffer() = delete;
  LogTreeBuffer(int log2size)
      : items(1 << log2size),
        present(1 << log2size, false),
        cur_size(0),
        insert_size(1 << log2size) {}

  size_t size() const { return cur_size; }
  bool empty() const { return cur_size == 0; }
  void clear() {
    insert_size = items.size();
    cur_size = 0;
    present.assign(items.size(), false);
  }

  size_t moveElementsTo(parlay::slice<objT *, objT *> dest) {
    auto cur_end = items.size() - insert_size;
    size_t ret;
    assert(dest.size() >= size());
    if (parallel) {
      ret = parlay::pack_into(parlay::slice(items.begin(), items.begin() + cur_end),
                              parlay::slice(present.begin(), present.begin() + cur_end),
                              dest);
    } else {
      ret = 0;
      for (size_t i = 0; i < cur_end; i++) {
        if (present[i]) {
          dest[ret++] = items[i];
        }
      }
    }
    clear();
    return ret;
  }

  void insert(const parlay::slice<const objT *, const objT *> &points) {
    if (points.size() == 0) return;
#ifndef NDEBUG
#warning "this is probably bad for runtime - remove later"
    if (points.size() + cur_size > items.size())
      throw std::runtime_error("invalid insertion into logtree buffer!");
#endif

    if (points.size() <= insert_size) {
      // place items
      auto cur_start = items.size() - insert_size;
      parlay::parallel_for(0, points.size(), [&](size_t i) {
        items[cur_start + i] = points[i];
        present[cur_start + i] = true;
      });
      // update size fields
      this->cur_size += points.size();
      this->insert_size -= points.size();
    } else {
      // gather elements
      parlay::sequence<objT> gather(cur_size);
      moveElementsTo(gather.cut(0, cur_size));
      // place elements
      auto new_size = gather.size() + points.size();
      parlay::parallel_for(0, new_size, [&](size_t i) {
        if (i < gather.size())
          this->items[i] = gather[i];
        else
          this->items[i] = points[i - gather.size()];
      });
      // update size fields
      this->cur_size = new_size;
      insert_size -= new_size;
    }
  }

  void bulk_erase(const parlay::sequence<objT> &points) {
    parlay::parallel_for(0, points.size(), [&](size_t i) {
      parlay::parallel_for(0, items.size() - insert_size, [&](size_t j) {
        if (present[j] && (points[i] == items[j])) {
          present[j] = false;
        }
      });
    });
  }

  template <bool _unused>
  void erase(const parlay::sequence<objT> &points) {
    bulk_erase(points);
  }

  // queries
  bool contains(const objT &p) const {
    // TODO: not great, has a data race
    bool ret = false;
    parlay::parallel_for(0, items.size() - insert_size, [&](size_t i) {
      if (present[i] && (p == items[i])) ret = true;
    });
    return ret;
  }

  parlay::sequence<objT> orthogonalQuery(const objT &qMin, const objT &qMax) const {
    auto cur_end = items.size() - insert_size;

    // check all points in parallel
    parlay::sequence<bool> in_box(cur_end, false);
    parlay::parallel_for(0, cur_end, [&](size_t i) {
      if (present[i]) {
        auto cmp = boxCompare(qMin, qMax, items[i], items[i]);
        if (cmp == BOX_INCLUDE) {
          in_box[i] = true;
        } else {
          assert(cmp == BOX_EXCLUDE);
        }
      }
    });

    // construct return sequence
    parlay::sequence<objT> ret(cur_end);
    auto ret_size =
        parlay::pack_into(parlay::slice(items.begin(), items.begin() + cur_end), in_box, ret);
    ret.resize(ret_size);
    return ret;
  }

  void knnSinglePoint(const objT &p, knnBuf::buffer<const pointT *> &buf) const {
    auto cur_end = items.size() - insert_size;
    for (size_t i = 0; i < cur_end; i++) {
      if (present[i]) {
        auto dist = p.dist(items[i]);
        // if (dist <= radius) {
        const objT *item_ptr = items.begin() + i;
        buf.insert(knnBuf::elem(dist, item_ptr));
        //}
      }
    }
  }

  template <bool set_res, bool _update, bool _recurse_sibling>
  void knnSinglePoint(
      const objT &p,
      int i,
      parlay::slice<knnBuf::elem<const pointT *> *, knnBuf::elem<const pointT *> *> &out,
      parlay::slice<const pointT **, const pointT **> &res,
      int k,
      bool preload) const {
    auto buf = knnBuf::buffer<const pointT *>(k, out.cut(i * 2 * k, (i + 1) * 2 * k));
    if (preload) buf.ptr = k;

    knnSinglePoint(p, buf);

    if (set_res) {
      for (int j = 0; j < k; j++) {
        res[i * k + j] = buf[j].entry;
      }
    }
  }

  void knn(const parlay::sequence<objT> &queries,
           parlay::slice<knnBuf::buffer<const pointT *> *, knnBuf::buffer<const pointT *> *> &bufs)
      const {
    if (parallel) {
      parlay::parallel_for(
          0, queries.size(), [&](size_t i) { knnSinglePoint(queries[i], bufs[i]); });
    } else {
      for (size_t i = 0; i < queries.size(); i++)
        knnSinglePoint(queries[i], bufs[i]);
    }
  }

  template <bool set_res, bool _update, bool _recurse_sibling>
  void knn(const parlay::sequence<objT> &queries,
           parlay::slice<knnBuf::elem<const pointT *> *, knnBuf::elem<const pointT *> *> &out,
           parlay::slice<const pointT **, const pointT **> &res,
           int k,
           bool preload = false) const {
    assert(res.size() == k * queries.size());
    assert(out.size() == 2 * k * queries.size());

    if (parallel) {
      parlay::parallel_for(0, queries.size(), [&](size_t i) {
        knnSinglePoint<set_res, false, false>(queries[i], i, out, res, k, preload);
      });
    } else {
      for (size_t i = 0; i < queries.size(); i++)
        knnSinglePoint<set_res, false, false>(queries[i], i, out, res, k, preload);
    }
  }
};

#endif  // LOGTREE_BUFFER_H
