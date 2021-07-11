#include "parlay/sequence.h"
#include "pargeo/point.h"

#include "convexHull2d/quickHull/hull.h"
#include "convexHull2d/divideConquer/hull.h"

template <class pointT>
inline typename pointT::floatT triArea(pointT a, pointT b, pointT c) {
  auto cross = [&](pointT p1, pointT p2) {
		 return p1[0] * p2[1] - p1[1] * p2[0];
	       };
  return cross((b-a), (c-a));
}

template<class pointT>
parlay::sequence<size_t>
pargeo::hull2d::divideConquer::compute(parlay::slice<pointT*, pointT*> P, size_t numProc) {
  if (!numProc) numProc = parlay::num_workers();

  numProc *= 8;

  size_t blkSize = floor(P.size() / numProc);

  while (blkSize < 10) {
    numProc -= 1;
    blkSize = floor(P.size() / numProc);
  }

  parlay::sequence<parlay::sequence<size_t>> subHulls(numProc);

  parlay::parallel_for(0, numProc, [&](size_t i) {
				     size_t s = i * blkSize;
				     size_t e = std::min(P.size(), (i+1) * blkSize);
				     //std::cout << s << " -- " << e << "\n";

				     parlay::sequence<size_t> h = std::move(pargeo::hull2d::quickHull::compute(P.cut(s, e)));
				     subHulls[i] = std::move(parlay::tabulate(h.size(),
									      [&](size_t i){
										return s + h[i];
									      }));

				   }, 1);

  parlay::sequence<size_t> uniquePts = parlay::flatten(subHulls);

  // point class with original index
  class pt : public pointT {
  public:
    pt(pointT p): pointT(p) {}
    pt(pointT p, typename pointT::attT att): pointT(p.coords(), att) {}
    size_t i;
  };

  parlay::sequence<pt> Q2 = parlay::tabulate(uniquePts.size(),
					     [&](size_t i){
					       auto tmp = pt(P[uniquePts[i]]);
					       tmp.i = uniquePts[i];
					       return tmp;
					     });

  auto H2 = pargeo::hull2d::quickHull::compute<pt>(parlay::make_slice(Q2));

  // restore original index
  return parlay::tabulate(H2.size(), [&](size_t i){
				       return Q2[H2[i]].i;
				     });
}

template parlay::sequence<size_t>
pargeo::hull2d::divideConquer::compute(parlay::slice<pargeo::fpoint<2>*, pargeo::fpoint<2>*>, size_t);

template parlay::sequence<size_t>
pargeo::hull2d::divideConquer::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*>, size_t);
