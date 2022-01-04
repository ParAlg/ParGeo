#pragma once

#include <math.h>
#include "pargeo/point.h"
#include "convexHull3d/facet.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/hash_table.h"

namespace pargeo {
  namespace hull3d {
    namespace samplingHelper {

      template<class pointT>
      parlay::sequence<pointT>
      randomSample(parlay::slice<pointT*, pointT*> P, size_t m) {
	return std::move(parlay::tabulate(m, [&](size_t i) {
					       return P[parlay::hash64(i) % P.size()];
					     }));
      }

      template<class pointT>
      parlay::sequence<pointT>
      randomProjection(parlay::slice<pointT*, pointT*> P, size_t m) {
	using floatT = typename pointT::floatT;
	m = m / 2;
	m = std::max(size_t(2), m);
	parlay::sequence<pointT> S(m*2);

	auto projection = [&](size_t i) {
			    size_t ii = parlay::hash64(i) % P.size();
			    pointT a = P[ii];
			    pointT b = P[P.size() - ii];
			    pointT ab = b - a;
			    auto minMax =
			      parlay::minmax_element(P.cut(ii, std::min(ii+m*2, P.size())),
						     [&](pointT p1, pointT p2) {
						       auto ap1 = p1 - a;
						       auto ap2 = p2 - a;
						       return (ap1.dot(ab) / ab.dot(ab)) <
							 (ap2.dot(ab) / ab.dot(ab));
						     });
			    S[i * 2] = minMax.first;
			    S[i * 2 + 1] = minMax.second;
			  };

	parlay::parallel_for(0, m, projection);
	return std::move(S);
      }

      template<class pointT>
      parlay::sequence<pointT>
      gridSample(parlay::slice<pointT*, pointT*> P, size_t m) {
	using namespace pargeo;
	using namespace parlay;
	using floatT = typename pointT::floatT;
	floatT d = pow(floatT(m), floatT(0.33333));
	size_t di = size_t(ceil(d));

	if (di > 2000000) // blows up size_t
	  throw std::runtime_error("too many grid samples");

	std::atomic<floatT> pMin[3];
	std::atomic<floatT> pMax[3];
	for (int i = 0; i < 3; ++ i) pMin[i] = std::numeric_limits<floatT>::max();
	for (int i = 0; i < 3; ++ i) pMax[i] = std::numeric_limits<floatT>::lowest();
	parallel_for(0, P.size(), [&](size_t i){
				    write_max(&pMax[0], P[i][0], std::less<floatT>());
				    write_min(&pMin[0], P[i][0], std::less<floatT>());
				    write_max(&pMax[1], P[i][1], std::less<floatT>());
				    write_min(&pMin[1], P[i][1], std::less<floatT>());
				    write_max(&pMax[2], P[i][2], std::less<floatT>());
				    write_min(&pMin[2], P[i][2], std::less<floatT>());
				  });

	floatT maxSpan = 0;
	for (int i = 0; i < 3; ++ i) maxSpan = std::max(maxSpan, pMax[i]-pMin[i]);
	floatT r = maxSpan / floatT(di);

	auto id = [&](pointT p) {
		    size_t x = floor((p[0] - pMin[0]) / r);
		    size_t y = floor((p[1] - pMin[1]) / r);
		    size_t z = floor((p[2] - pMin[2]) / r);
		    return x*di*di + y*di + z;
		  };

	hashtable<hash_numeric<size_t>> T(P.size(), hash_numeric<size_t>());
	auto sample = tabulate(P.size(), [&](size_t i){
					   if (T.insert(id(P[i]))) return P[i];
					   else return pointT();
					 });

	return parlay::filter(sample, [&](pointT p) {return !p.isEmpty();});
      }

    } // End namespace sampling
  } // End namespace hullInternal
} // End namespace pargeo
