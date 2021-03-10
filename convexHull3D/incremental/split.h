#pragma once

#include <iostream>

#include "parlay/delayed_sequence.h"
#include "parlay/monoid.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"

using namespace parlay;
using namespace parlay::internal;

template <typename InIterator, typename OutIterator, typename Char_Seq>
sequence<size_t> split_k(size_t k, slice<InIterator, InIterator> In,
			 slice<OutIterator, OutIterator> Out,
			 Char_Seq const &Fl,
			 flags fl = no_flag) {
  size_t n = In.size();
  if (In == Out)
    throw std::invalid_argument("In and Out cannot be the same in split_four");
  size_t l = num_blocks(n, _block_size);

  sequence<size_t>* Sums[k-1];
  for (size_t i=0; i<k-1; ++i) Sums[i] = new sequence<size_t>(l);

  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
	       size_t c[k-1];
	       for (size_t x=0; x<k-1; ++x) c[x] = 0;

               for (size_t j = s; j < e; j++) {
		 for (size_t x=0; x<k-1; ++x) {
		   if (Fl[j] == x) c[x]++;
		 }
               }
	       for (size_t x=0; x<k-1; ++x)
		 Sums[x]->at(i) = c[x];
             },
             fl);
  size_t m[k-1];
  for (size_t x=0; x<k-1; ++x)
    m[x] = scan_inplace(Sums[x]->cut(0, Sums[x]->size()), addm<size_t>());
  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
	       size_t c[k];
	       c[0] = Sums[0]->at(i);
	       c[1] = m[0] + Sums[1]->at(i);
	       for (size_t x=2; x<k; ++x) {
		 c[x] = s;
		 for (size_t y=0; y<x; ++y) {
		   c[x] += m[y];
		   c[x] -= Sums[y]->at(i);}
	       }
               for (size_t j = s; j < e; j++) {
		 for (size_t x=0; x<k; x++) {
		   if (Fl[j] == x) Out[c[x]++] = In[j];
		 }
               }
             },
             fl);
  for (size_t x=0; x<k-1; x++) delete Sums[x];
  auto M = sequence<size_t>(k-1);
  for (size_t x=0; x<k-1; x++)
    M[x] = m[x];
  return M;
}
