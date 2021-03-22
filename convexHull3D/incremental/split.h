#pragma once

#include <iostream>

#include "parlay/delayed_sequence.h"
#include "parlay/monoid.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"

using namespace parlay;
using namespace parlay::internal;

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
