#pragma once

#include <math.h>
#include <atomic>
#include "pargeo/getTime.h"
#include "kdTree/kdTree.h"
#include "pargeo/point.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/hash_table.h"

// #define SAMPLE_HULL_VERBOSE

namespace pargeo {
  namespace hull2d {
    namespace samplingHelper {

      template <class pt>
      class facet {
	using floatT = typename pt::floatT;

      private:

	inline floatT triArea(pt a, pt b, pt c) {
	  auto cross = [&](pt p1, pt p2) {
	    return p1[0] * p2[1] - p1[1] * p2[0];
	  };
	  return cross((b-a), (c-a));
	}

      public:
	pt a, b;

	facet(pt _a, pt _b): a(_a), b(_b) {}

	floatT getArea(pt c) {
	  return triArea(a, b, c);
	}

      };

      template<class pointT>
      std::tuple<typename pointT::floatT, typename pointT::floatT>
      nodeVol(kdTree::node<2, pointT>* r, facet<pointT>* f) {
	using floatT = typename pointT::floatT;

	floatT vMin = std::numeric_limits<floatT>::max();
	floatT vMax = -1;
	floatT vol;
	pointT c;

	c[0] = r->getMin(0); c[1] = r->getMin(1);
	vol = f->getArea(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);

	c[0] = r->getMin(0); c[1] = r->getMax(1);
	vol = f->getArea(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);

	c[0] = r->getMax(0); c[1] = r->getMin(1);
	vol = f->getArea(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);

	c[0] = r->getMax(0); c[1] = r->getMax(1);
	vol = f->getArea(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);

	return std::tuple(vMin, vMax);
      }

      template<class pointT>
      void flagVisible(kdTree::node<2, pointT>* r,
		       facet<pointT>* f,
		       parlay::slice<std::atomic<bool>*, std::atomic<bool>*> flag,
		       pointT* base) {
	using floatT = typename pointT::floatT;
	floatT eps = pointT::eps;

	auto vols = nodeVol<pointT>(r, f);
	floatT volMin = std::get<0>(vols);
	floatT volMax = std::get<1>(vols);
	bool F = false;

	if (volMin <= 0 && volMax <= 0) {
	  // node is in hull
	  return;
	} else if (volMin <= 0 && volMax > 0) {
	  // node intersect
	  if (!r->isLeaf()) {
	    flagVisible<pointT>(r->L(), f, flag, base);
	    flagVisible<pointT>(r->R(), f, flag, base);
	  } else {
	    for (size_t i = 0; i < r->size(); ++ i) {
	      auto pVol = f->getArea(*r->at(i));
	      if (pVol >= 0 ||
		  f->a == *r->at(i) ||
		  f->b == *r->at(i)
		  ) {
		std::atomic_compare_exchange_weak(&flag[r->at(i) - base], &F, true);
	      }
	    }
	  }
	} else {
	  // node is outside of hull
	  for (size_t i = 0; i < r->size(); ++ i) {
	    std::atomic_compare_exchange_weak(&flag[r->at(i) - base], &F, true);
	  }
	}

      }

      template<class pointT>
      parlay::sequence<pointT>
      filter2(parlay::slice<pointT*, pointT*> P,
	      parlay::slice<pointT*, pointT*> sample,
	      parlay::slice<size_t*, size_t*> F) {
	using namespace parlay;

	/* turn hullIndices F into facets */

	sequence<facet<pointT>> facets = parlay::tabulate(F.size(), [&](size_t i) {
	  return facet(sample[F[i]], sample[F[(i+1) % F.size()]]);
	});

	/* checking code, just filter here directly without the tree */
	// parlay::sequence<pointT> Q;
	// for (auto p: P) {
	//   for (auto f: facets) {
	//     if (f.getArea(p) >= p.eps || p == f.a || p == f.b) {
	//       Q.push_back(p);
	//       break;
	//     }}
	// }
	// return Q;

#ifdef SAMPLE_HULL_VERBOSE
	pargeo::timer t; t.start();
#endif

	pargeo::kdTree::node<2, pointT>* tree =
	  kdTree::build<2, pointT>(P, true, 128); // todo manual leaf size

#ifdef SAMPLE_HULL_VERBOSE
	std::cout << " build-tree-time = " << t.get_next() << "\n";
#endif

	parlay::sequence<std::atomic<bool>> flag(P.size());
	parlay::parallel_for(0, flag.size(), [&](size_t i){flag[i] = false;});

	auto interiorPt = (facets[0].a + facets[facets.size()-1].a) / 2;

	auto data = &P[0];
	parlay::parallel_for(0, F.size(), [&](size_t i) {
	  flagVisible<pointT>(tree, &facets[i], make_slice(flag), data);
	});

	pargeo::kdTree::del(tree);

#ifdef SAMPLE_HULL_VERBOSE
	std::cout << " filter-time = " << t.get_next() << "\n";
#endif

	auto POut = std::move(parlay::pack(P, flag));

#ifdef SAMPLE_HULL_VERBOSE
	std::cout << " pack-time = " << t.get_next() << "\n";
#endif
	return std::move(POut);
      }

    } // End namespace helper
  } // End namespace hullInternal
} // End namespace pargeo
