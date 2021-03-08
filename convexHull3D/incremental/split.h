#ifndef SPLIT_NINE_H
#define SPLIT_NINE_H

#include <iostream>

#include "parlay/delayed_sequence.h"
#include "parlay/monoid.h"
#include "parlay/sequence.h"
#include "parlay/utilities.h"

using namespace parlay;
using namespace parlay::internal;

template <typename InIterator, typename OutIterator, typename Char_Seq>
sequence<size_t> split_four(slice<InIterator, InIterator> In,
                                      slice<OutIterator, OutIterator> Out,
                                      Char_Seq const &Fl, flags fl = no_flag) {
  size_t n = In.size();
  if (In == Out)
    throw std::invalid_argument("In and Out cannot be the same in split_four");
  size_t l = num_blocks(n, _block_size);
  sequence<size_t> Sums0(l);
  sequence<size_t> Sums1(l);
  sequence<size_t> Sums2(l);
  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c0 = 0;
               size_t c1 = 0;
               size_t c2 = 0;
               for (size_t j = s; j < e; j++) {
                 if (Fl[j] == 0) c0++;
                 else if (Fl[j] == 1) c1++;
                 else if (Fl[j] == 2) c2++;
               }
               Sums0[i] = c0;
               Sums1[i] = c1;
               Sums2[i] = c2;
             },
             fl);
  size_t m0 = scan_inplace(make_slice(Sums0), addm<size_t>());
  size_t m1 = scan_inplace(make_slice(Sums1), addm<size_t>());
  size_t m2 = scan_inplace(make_slice(Sums2), addm<size_t>());
  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c0 = Sums0[i];
               size_t c1 = m0 + Sums1[i];
               size_t c2 = m0 + m1 + (s - Sums0[i] - Sums1[i]);
               size_t c3 = m0 + m1 + m2 + (s - Sums0[i] - Sums1[i] - Sums2[i]);
               for (size_t j = s; j < e; j++) {
                 if (Fl[j] == 0) Out[c0++] = In[j];
                 else if (Fl[j] == 1) Out[c1++] = In[j];
                 else if (Fl[j] == 2) Out[c2++] = In[j];
                 else Out[c3++] = In[j];
               }
             },
             fl);
  auto M = sequence<size_t>(3);
  M[0] = m0;
  M[1] = m1;
  M[2] = m2;
  return M;
}

template <typename InIterator, typename OutIterator, typename Char_Seq>
sequence<size_t> split_nine(slice<InIterator, InIterator> In,
                                      slice<OutIterator, OutIterator> Out,
                                      Char_Seq const &Fl, flags fl = no_flag) {
  size_t n = In.size();
  if (In == Out)
    throw std::invalid_argument("In and Out cannot be the same in split_nine");
  size_t l = num_blocks(n, _block_size);
  sequence<size_t> Sums0(l);
  sequence<size_t> Sums1(l);
  sequence<size_t> Sums2(l);
  sequence<size_t> Sums3(l);
  sequence<size_t> Sums4(l);
  sequence<size_t> Sums5(l);
  sequence<size_t> Sums6(l);
  sequence<size_t> Sums7(l);
  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c0 = 0;
               size_t c1 = 0;
               size_t c2 = 0;
               size_t c3 = 0;
               size_t c4 = 0;
               size_t c5 = 0;
               size_t c6 = 0;
               size_t c7 = 0;
               for (size_t j = s; j < e; j++) {
                 if (Fl[j] == 0) c0++;
                 else if (Fl[j] == 1) c1++;
                 else if (Fl[j] == 2) c2++;
                 else if (Fl[j] == 3) c3++;
                 else if (Fl[j] == 4) c4++;
                 else if (Fl[j] == 5) c5++;
                 else if (Fl[j] == 6) c6++;
                 else if (Fl[j] == 7) c7++;
               }
               Sums0[i] = c0;
               Sums1[i] = c1;
               Sums2[i] = c2;
               Sums3[i] = c3;
               Sums4[i] = c4;
               Sums5[i] = c5;
               Sums6[i] = c6;
               Sums7[i] = c7;
             },
             fl);
  size_t m0 = scan_inplace(make_slice(Sums0), addm<size_t>());
  size_t m1 = scan_inplace(make_slice(Sums1), addm<size_t>());
  size_t m2 = scan_inplace(make_slice(Sums2), addm<size_t>());
  size_t m3 = scan_inplace(make_slice(Sums3), addm<size_t>());
  size_t m4 = scan_inplace(make_slice(Sums4), addm<size_t>());
  size_t m5 = scan_inplace(make_slice(Sums5), addm<size_t>());
  size_t m6 = scan_inplace(make_slice(Sums6), addm<size_t>());
  size_t m7 = scan_inplace(make_slice(Sums7), addm<size_t>());
  sliced_for(n, _block_size,
             [&](size_t i, size_t s, size_t e) {
               size_t c0 = Sums0[i];
               size_t c1 = m0 + Sums1[i];
               size_t c2 = m0 + m1 + (s - Sums0[i] - Sums1[i]);
               size_t c3 = m0 + m1 + m2 + (s - Sums0[i] - Sums1[i] - Sums2[i]);
               size_t c4 = m0 + m1 + m2 + m3 + (s - Sums0[i] - Sums1[i] - Sums2[i] - Sums3[i]);
               size_t c5 = m0 + m1 + m2 + m3 + m4 + (s - Sums0[i] - Sums1[i] - Sums2[i] - Sums3[i] - Sums4[i]);
               size_t c6 = m0 + m1 + m2 + m3 + m4 + m5 + (s - Sums0[i] - Sums1[i] - Sums2[i] - Sums3[i] - Sums4[i] - Sums5[i]);
               size_t c7 = m0 + m1 + m2 + m3 + m4 + m5 + m6 + (s - Sums0[i] - Sums1[i] - Sums2[i] - Sums3[i] - Sums4[i] - Sums5[i] - Sums6[i]);
               size_t c8 = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + (s - Sums0[i] - Sums1[i] - Sums2[i] - Sums3[i] - Sums4[i] - Sums5[i] - Sums6[i] - Sums7[i]);
               for (size_t j = s; j < e; j++) {
                 if (Fl[j] == 0) Out[c0++] = In[j];
                 else if (Fl[j] == 1) Out[c1++] = In[j];
                 else if (Fl[j] == 2) Out[c2++] = In[j];
                 else if (Fl[j] == 3) Out[c3++] = In[j];
                 else if (Fl[j] == 4) Out[c4++] = In[j];
                 else if (Fl[j] == 5) Out[c5++] = In[j];
                 else if (Fl[j] == 6) Out[c6++] = In[j];
                 else if (Fl[j] == 7) Out[c7++] = In[j];
                 else Out[c8++] = In[j];
               }
             },
             fl);
  auto M = sequence<size_t>(8);
  M[0] = m0;
  M[1] = m1;
  M[2] = m2;
  M[3] = m3;
  M[4] = m4;
  M[5] = m5;
  M[6] = m6;
  M[7] = m7;
  return M;
}

#endif
