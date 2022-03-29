// This code is part of the project "Pargeo: A Library for Parallel Computational Geometry"
// Copyright (c) 2022 Yiqiu Wang, Shangdi Yu, Laxman Dhulipala, Yan Gu, Julian Shun
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

#include "mortonSort/mortonSort.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/point.h"
#include <tuple>

// This version is customized for batchKdTree, and it only works with
// the point class in batchKdTree/shared/geometry.h
template<typename pt>
parlay::sequence<pt> pargeo::zorderSort2d_2(const parlay::sequence<pt>& P) {

  using floatT = typename pt::floatT;

  static constexpr size_t maxRange = 65535;
  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  auto shift = [&](size_t x) {
		 if (x > maxRange) {
		   //cout << x << endl;
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
		extrema[i*2] = P[0].readCoord(i);
		extrema[i*2+1] = P[0].readCoord(i);
  }

  parlay::parallel_for(0, P.size(), [&](size_t i){
			      parlay::write_max(&extrema[0], P[i].readCoord(0), std::less<floatT>());
			      parlay::write_min(&extrema[1], P[i].readCoord(0), std::less<floatT>());
			      parlay::write_max(&extrema[2], P[i].readCoord(1), std::less<floatT>());
			      parlay::write_min(&extrema[3], P[i].readCoord(1), std::less<floatT>());
			    });

  floatT maxSpan = std::max((extrema[2]-extrema[3]),(extrema[0]-extrema[1]));

  floatT boxSize = 1.01 * maxSpan / maxRange;

  point<2> pMin;
  pMin[0] = extrema[1];
  pMin[1] = extrema[3];

  using ip = std::tuple<size_t, pt>;
  parlay::sequence<ip> pairs = parlay::tabulate(P.size(), [&](size_t i){
					    return std::tuple(zorder(P[i], pMin, boxSize), P[i]);
					  });

  parlay::sort_inplace(parlay::make_slice(pairs), [&](ip i, ip j){
						    return std::get<0>(i) < std::get<0>(j);
					  });

  return parlay::tabulate(pairs.size(), [&](size_t i){
				  return std::get<1>(pairs[i]);
				});

	return P;
}

template<typename pt>
parlay::sequence<pt> pargeo::zorderSort2d(parlay::sequence<pt>& P) {

  using floatT = typename pt::floatT;

  static constexpr size_t maxRange = 65535;
  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  auto shift = [&](size_t x) {
		 if (x > maxRange) {
		   //cout << x << endl;
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

  parlay::parallel_for(0, P.size(), [&](size_t i){
			      parlay::write_max(&extrema[0], P[i][0], std::less<floatT>());
			      parlay::write_min(&extrema[1], P[i][0], std::less<floatT>());
			      parlay::write_max(&extrema[2], P[i][1], std::less<floatT>());
			      parlay::write_min(&extrema[3], P[i][1], std::less<floatT>());
			    });

  floatT maxSpan = std::max((extrema[2]-extrema[3]),(extrema[0]-extrema[1]));

  floatT boxSize = 1.01 * maxSpan / maxRange;

  point<2> pMin;
  pMin[0] = extrema[1];
  pMin[1] = extrema[3];

  using ip = std::tuple<size_t, pt>;
  parlay::sequence<ip> pairs = parlay::tabulate(P.size(), [&](size_t i){
					    return std::tuple(zorder(P[i], pMin, boxSize), P[i]);
					  });

  parlay::sort_inplace(parlay::make_slice(pairs), [&](ip i, ip j){
						    return std::get<0>(i) < std::get<0>(j);
					  });

  return parlay::tabulate(pairs.size(), [&](size_t i){
				  return std::get<1>(pairs[i]);
				});
}

// This version is customized for batchKdTree, and it only works with
// the point class in batchKdTree/shared/geometry.h
template<typename pt>
parlay::sequence<pt> pargeo::zorderSort3d_2(const parlay::sequence<pt>& P) {

  using floatT = typename pt::floatT;

  static constexpr size_t maxRange = 65535;
  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  auto shift = [&](size_t x) {
		 if (x > maxRange) {
		   //cout << x << endl;
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
		extrema[i*2] = P[0].readCoord(i);
    extrema[i*2+1] = P[0].readCoord(i);
  }

  parlay::parallel_for(0, P.size(),
		       [&](size_t i){
			 parlay::write_max(&extrema[0], P[i].readCoord(0), std::less<floatT>());
			 parlay::write_min(&extrema[1], P[i].readCoord(0), std::less<floatT>());
			 parlay::write_max(&extrema[2], P[i].readCoord(1), std::less<floatT>());
			 parlay::write_min(&extrema[3], P[i].readCoord(1), std::less<floatT>());
			 parlay::write_max(&extrema[4], P[i].readCoord(2), std::less<floatT>());
			 parlay::write_min(&extrema[5], P[i].readCoord(2), std::less<floatT>());
	});

  floatT maxSpan = std::max(extrema[4]-extrema[5],
			    std::max((extrema[2]-extrema[3]),
				     (extrema[0]-extrema[1])));

  floatT boxSize = 1.01 * maxSpan / maxRange;

  point<3> pMin;
  pMin[0] = extrema[1];
  pMin[1] = extrema[3];
  pMin[2] = extrema[5];

  using ip = std::tuple<size_t, pt>;
  parlay::sequence<ip> pairs = parlay::tabulate(P.size(), [&](size_t i){
					    return std::tuple(zorder(P[i], pMin, boxSize), P[i]);
					  });

  parlay::sort_inplace(parlay::make_slice(pairs), [&](ip i, ip j){
						    return std::get<0>(i) < std::get<0>(j);
					  });

  return parlay::tabulate(pairs.size(), [&](size_t i){
				  return std::get<1>(pairs[i]);
				});
}

template<typename pt>
parlay::sequence<pt> pargeo::zorderSort3d(parlay::sequence<pt>& P) {

  using floatT = typename pt::floatT;

  static constexpr size_t maxRange = 65535;
  static constexpr size_t maxBit = 16; // 2**maxBit = maxRange

  auto shift = [&](size_t x) {
		 if (x > maxRange) {
		   //cout << x << endl;
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

  parlay::parallel_for(0, P.size(),
		       [&](size_t i){
			 parlay::write_max(&extrema[0], P[i][0], std::less<floatT>());
			 parlay::write_min(&extrema[1], P[i][0], std::less<floatT>());
			 parlay::write_max(&extrema[2], P[i][1], std::less<floatT>());
			 parlay::write_min(&extrema[3], P[i][1], std::less<floatT>());
			 parlay::write_max(&extrema[4], P[i][2], std::less<floatT>());
			 parlay::write_min(&extrema[5], P[i][2], std::less<floatT>());
		       });

  floatT maxSpan = std::max(extrema[4]-extrema[5],
			    std::max((extrema[2]-extrema[3]),
				     (extrema[0]-extrema[1])));

  floatT boxSize = 1.01 * maxSpan / maxRange;

  point<3> pMin;
  pMin[0] = extrema[1];
  pMin[1] = extrema[3];
  pMin[2] = extrema[5];

  using ip = std::tuple<size_t, pt>;
  parlay::sequence<ip> pairs = parlay::tabulate(P.size(), [&](size_t i){
					    return std::tuple(zorder(P[i], pMin, boxSize), P[i]);
					  });

  parlay::sort_inplace(parlay::make_slice(pairs), [&](ip i, ip j){
						    return std::get<0>(i) < std::get<0>(j);
					  });

  return parlay::tabulate(pairs.size(), [&](size_t i){
				  return std::get<1>(pairs[i]);
				});
}

template<typename pt>
void pargeo::zorderSortInPlace2d(parlay::sequence<pt>& P) {
  P = zorderSort2d(P);
}

template<typename pt>
void pargeo::zorderSortInPlace3d(parlay::sequence<pt>& P) {
  P = zorderSort3d(P);
}
