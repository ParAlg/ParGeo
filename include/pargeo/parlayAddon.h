#pragma once

#include "parlay/sequence.h"
#include "parlay/primitives.h"
#include "parlay/internal/sequence_ops.h"

namespace parlay {
  using namespace parlay::internal;

  /* -------------------- Serial Min and max -------------------- */

  template <PARLAY_RANGE_TYPE R, typename Compare>
    auto minmax_element_serial(R&& r, Compare comp) {
    size_t n = parlay::size(r);
    auto SS = delayed_seq<std::pair<size_t, size_t>>(parlay::size(r),
                [&](size_t i) { return std::make_pair(i, i); });
    auto f = [&comp, it = std::begin(r)](const auto& l, const auto& r) {
               return (std::make_pair(!comp(it[r.first], it[l.first]) ? l.first : r.first,
                 !comp(it[l.second], it[r.second]) ? l.second : r.second));
    };
    auto ds = internal::reduce(make_slice(SS), make_monoid(f, std::make_pair(n, n)));
    return std::make_pair(std::begin(r) + ds.first, std::begin(r) + ds.second);
  }

  template <PARLAY_RANGE_TYPE R, typename Compare>
    auto min_element_serial(R&& r, Compare comp) {
    auto SS = delayed_seq<size_t>(parlay::size(r), [&](size_t i) { return i; });
    auto f = [&comp, it = std::begin(r)](size_t l, size_t r)
      { return (!comp(it[r], it[l]) ? l : r); };
    return std::begin(r) +
      internal::reduce_serial(make_slice(SS), make_monoid(f, (size_t)parlay::size(r)));
  }

  template <PARLAY_RANGE_TYPE R>
    auto min_element_serial(R&& r) {
    return min_element_serial(r, std::less<range_value_type_t<R>>());
  }

  template <PARLAY_RANGE_TYPE R, typename Compare>
    auto max_element_serial(R&& r, Compare comp) {
    return min_element_serial(r, [&](const auto& a, const auto& b)
		       { return comp(b, a); });
  }

  template <PARLAY_RANGE_TYPE R>
    auto max_element_serial(R&& r) {
    return max_element_serial(r, std::less<range_value_type_t<R>>());
  }

  /* -------------------- K-way splitting functions -------------------- */

  template <typename mySeq, typename boolSeq>
  sequence<mySeq> split_k_2(size_t k, mySeq& In,
			   boolSeq const &Fl,
			   flags fl = no_flag) {
    size_t n = In.size();
    size_t l = num_blocks(n, _block_size);

    sequence<size_t>* Sums[k];
    for (size_t i=0; i<k; ++i) Sums[i] = new sequence<size_t>(l);

    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c[k];
		 for (size_t x=0; x<k; ++x) c[x] = 0;

		 for (size_t j = s; j < e; j++) {
		   for (size_t x=0; x<k; ++x) {
		     if (Fl[j] == x) c[x]++;
		   }
		 }
		 for (size_t x=0; x<k; ++x)
		   Sums[x]->at(i) = c[x];
	       },
	       fl);

    auto chunks = sequence<mySeq>(k);
    size_t m[k];
    for (size_t x=0; x<k; ++x) {
      m[x] = scan_inplace(Sums[x]->cut(0, Sums[x]->size()), addm<size_t>());
      chunks[x] = mySeq(m[x]);
    }

    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c[k];
		 for (size_t x = 0; x < k; ++x) c[x] = 0;

		 for (size_t j = s; j < e; j++) {
		   for (size_t x = 0; x < k; x++) {
		     if (Fl[j] == x) {
		       chunks[x][Sums[x]->at(i) + c[x]++] = In.at(j);
		       break;
		     }
		   }
		 }
	       },
	       fl);
    for (size_t x=0; x<k; x++) delete Sums[x];

    return chunks;
  }

  template <typename mySeq, typename boolSeq>
  sequence<mySeq*> split_k(size_t k, mySeq* In,
			   boolSeq const &Fl,
			   flags fl = no_flag) {
    size_t n = In->size();
    size_t l = num_blocks(n, _block_size);

    sequence<size_t>* Sums[k];
    for (size_t i=0; i<k; ++i) Sums[i] = new sequence<size_t>(l);

    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c[k];
		 for (size_t x=0; x<k; ++x) c[x] = 0;

		 for (size_t j = s; j < e; j++) {
		   for (size_t x=0; x<k; ++x) {
		     if (Fl[j] == x) c[x]++;
		   }
		 }
		 for (size_t x=0; x<k; ++x)
		   Sums[x]->at(i) = c[x];
	       },
	       fl);

    auto chunks = sequence<mySeq*>(k);
    size_t m[k];
    for (size_t x=0; x<k; ++x) {
      m[x] = scan_inplace(Sums[x]->cut(0, Sums[x]->size()), addm<size_t>());
      chunks[x] = new mySeq(m[x]);
    }

    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c[k];
		 for (size_t x = 0; x < k; ++x) c[x] = 0;

		 for (size_t j = s; j < e; j++) {
		   for (size_t x = 0; x < k; x++) {
		     if (Fl[j] == x) {
		       chunks[x]->at(Sums[x]->at(i) + c[x]++) = In->at(j);
		       break;
		     }
		   }
		 }
	       },
	       fl);
    for (size_t x=0; x<k; x++) delete Sums[x];

    return chunks;
  }

} // End namespace
