// This code is part of the "Pargeo" project
// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "point.h"
#include <tuple>

namespace pargeo {

  using namespace parlay;
  using namespace std;

  template<typename pt>
  parlay::sequence<pt> zorderSort2d(parlay::sequence<pt>& P) {

    using floatT = double;

    static constexpr size_t maxRange = 65535;
    static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

    auto shift = [&](size_t x) {
		   if (x > maxRange) {
		     cout << x << endl;
		     throw std::runtime_error("zorder range error");
		   }

		   x = (x | (x << 16)) & 0x001f0000ff0000ff;
		   x = (x | (x <<  8)) & 0x100f00f00f00f00f;
		   x = (x | (x <<  4)) & 0x10c30c30c30c30c3;
		   x = (x | (x <<  2)) & 0x1249249249249249;
		   return x;
		 };

    auto zorder = [&](pt p, point<2> pMin, floatT boxSize) {
		    size_t x = floor( (p[0] - pMin[0]) / boxSize );
		    size_t y = floor( (p[1] - pMin[1]) / boxSize );

		    x = shift(x);
		    y = shift(y);
		    return x | (y << 1);
		  };

    std::atomic<floatT> extrema[4];

    for (int i=0; i<2; ++i) {
      extrema[i*2] = P[0][i];
      extrema[i*2+1] = P[0][i];
    }

    parallel_for(0, P.size(), [&](size_t i){
				 write_max(&extrema[0], P[i][0], std::less<floatT>());
				 write_min(&extrema[1], P[i][0], std::less<floatT>());
				 write_max(&extrema[2], P[i][1], std::less<floatT>());
				 write_min(&extrema[3], P[i][1], std::less<floatT>());
			       });

    floatT maxSpan = max((extrema[2]-extrema[3]),(extrema[0]-extrema[1]));

    floatT boxSize = 1.01 * maxSpan / maxRange;

    point<2> pMin;
    pMin[0] = extrema[1];
    pMin[1] = extrema[3];

    using ip = tuple<size_t, pt>;
    sequence<ip> pairs = tabulate(P.size(), [&](size_t i){
			   return tuple(zorder(P[i], pMin, boxSize), P[i]);
			 });

    parlay::sort_inplace(make_slice(pairs), [&](ip i, ip j){
					    return get<0>(i) < get<0>(j);
					  });

    return tabulate(pairs.size(), [&](size_t i){
				    return get<1>(pairs[i]);
				  });
  }

  template<typename pt>
  parlay::sequence<pt> zorderSort3d(parlay::sequence<pt>& P) {

    using floatT = double;

    static constexpr size_t maxRange = 65535;
    static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

    auto shift = [&](size_t x) {
		   if (x > maxRange) {
		     cout << x << endl;
		     throw std::runtime_error("zorder range error");
		   }

		   x = (x | (x << 16)) & 0x001f0000ff0000ff;
		   x = (x | (x <<  8)) & 0x100f00f00f00f00f;
		   x = (x | (x <<  4)) & 0x10c30c30c30c30c3;
		   x = (x | (x <<  2)) & 0x1249249249249249;
		   return x;
		 };

    auto zorder = [&](pt p, point<3> pMin, floatT boxSize) {
		    size_t x = floor( (p[0] - pMin[0]) / boxSize );
		    size_t y = floor( (p[1] - pMin[1]) / boxSize );
		    size_t z = floor( (p[2] - pMin[2]) / boxSize );

		    x = shift(x);
		    y = shift(y);
		    z = shift(z);
		    return x | (y << 1) | (z << 2);
		  };

    std::atomic<floatT> extrema[6];

    for (int i=0; i<3; ++i) {
      extrema[i*2] = P[0][i];
      extrema[i*2+1] = P[0][i];
    }

    parallel_for(0, P.size(), [&](size_t i){
				 write_max(&extrema[0], P[i][0], std::less<floatT>());
				 write_min(&extrema[1], P[i][0], std::less<floatT>());
				 write_max(&extrema[2], P[i][1], std::less<floatT>());
				 write_min(&extrema[3], P[i][1], std::less<floatT>());
				 write_max(&extrema[4], P[i][2], std::less<floatT>());
				 write_min(&extrema[5], P[i][2], std::less<floatT>());
			       });

    floatT maxSpan = max(extrema[4]-extrema[5],
		  max((extrema[2]-extrema[3]),(extrema[0]-extrema[1])));

    floatT boxSize = 1.01 * maxSpan / maxRange;

    point<3> pMin;
    pMin[0] = extrema[1];
    pMin[1] = extrema[3];
    pMin[2] = extrema[5];

    using ip = tuple<size_t, pt>;
    sequence<ip> pairs = tabulate(P.size(), [&](size_t i){
			   return tuple(zorder(P[i], pMin, boxSize), P[i]);
			 });

    parlay::sort_inplace(make_slice(pairs), [&](ip i, ip j){
					    return get<0>(i) < get<0>(j);
					  });

    return tabulate(pairs.size(), [&](size_t i){
				    return get<1>(pairs[i]);
				  });
  }

} // End namespace
