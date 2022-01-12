#pragma once

#include <math.h>
#include <atomic>
#include "pargeo/getTime.h"
#include "kdTree/kdTree.h"
#include "pargeo/point.h"
#include "convexHull3d/facet.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/hash_table.h"

// #define SAMPLE_HULL_VERBOSE

namespace pargeo {
  namespace hull3d {
    namespace samplingHelper {

      template<class pointT>
      parlay::sequence<pointT>
      filter1(parlay::slice<pointT*, pointT*> P,
	      parlay::slice<
	      pargeo::hull3d::facet<pointT>*,
	      pargeo::hull3d::facet<pointT>*> F) {

	auto interiorPt = (F[0].a + F[F.size()-1].a) / 2;

	parlay::sequence<pointT> area(F.size());
	parlay::parallel_for(0, F.size(),
		     [&](size_t i){
		       auto f = F[i];
		       area[i] = pargeo::crossProduct3d(f.b-f.a, f.c-f.a);
		     });

	auto isOut = [&](pointT p) {
		       for (size_t i = 0; i < F.size(); ++i) {
			 auto f = F[i];
			 if (((f.a - interiorPt) - (p - interiorPt)).dot(area[i]) > 0
			     || f.a == p || f.b == p || f.c == p
			     ) {
			   return true;
			 }
		       }
		       return false;
		     };

	return parlay::filter(parlay::make_slice(P), isOut);
      }

      template<class pt>
      struct facet : public pargeo::hull3d::facet<pt> {
      public:

	using floatT = typename pt::floatT;
	using T = pargeo::hull3d::facet<pt>;

	pt area;

	floatT getVolume(pt p) {
	  return (T::a - p).dot(area);
	}

	facet(pt _a, pt _b, pt _c): T(_a, _b, _c) {
	  area = pargeo::crossProduct3d(T::b - T::a, T::c - T::a);
	}

      };

      template<class pointT>
      std::tuple<typename pointT::floatT, typename pointT::floatT>
      nodeVol(kdNode<3, pointT>* r, facet<pointT>* f) {
	using floatT = typename pointT::floatT;

	floatT vMin = std::numeric_limits<floatT>::max();
	floatT vMax = -1;
	floatT vol;
	pointT c;

	c[0] = r->getMin(0); c[1] = r->getMin(1); c[2] = r->getMin(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMax(0); c[1] = r->getMin(1); c[2] = r->getMin(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMin(0); c[1] = r->getMax(1); c[2] = r->getMin(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMin(0); c[1] = r->getMin(1); c[2] = r->getMax(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMax(0); c[1] = r->getMax(1); c[2] = r->getMin(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMin(0); c[1] = r->getMax(1); c[2] = r->getMax(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMax(0); c[1] = r->getMin(1); c[2] = r->getMax(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	c[0] = r->getMax(0); c[1] = r->getMax(1); c[2] = r->getMax(2);
	vol = f->getVolume(c);
	vMax = std::max(floatT(vMax), vol); vMin = std::min(floatT(vMin), vol);
	return std::tuple(vMin, vMax);
      }

      template<class pointT>
      void flagVisible(kdNode<3, pointT>* r,
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
	      auto pVol = f->getVolume(*r->at(i));
	      if (pVol >= 0 ||
		  f->a == *r->at(i) ||
		  f->b == *r->at(i) ||
		  f->c == *r->at(i)
		  ) {
		//flag[r->at(i) - base] = true;
		std::atomic_compare_exchange_weak(&flag[r->at(i) - base], &F, true);
	      }
	    }
	  }
	} else {
	  // node is outside of hull
	  for (size_t i = 0; i < r->size(); ++ i) {
	    //auto pVol = f->getVolume(*r->at(i));
	    // if (pVol >= 0 ||
	    // 	f->a == *r->at(i) ||
	    // 	f->b == *r->at(i) ||
	    // 	f->c == *r->at(i)
	    // 	) {
	    std::atomic_compare_exchange_weak(&flag[r->at(i) - base], &F, true);
	    //flag[r->at(i) - base] = true;
	    //}
	  }
	}

      }

      template<class pointT>
      parlay::sequence<pointT>
      filter2(parlay::slice<pointT*, pointT*> P,
	      parlay::slice<
	      pargeo::hull3d::facet<pointT>*,
	      pargeo::hull3d::facet<pointT>*> F) {
	using namespace parlay;

#ifdef SAMPLE_HULL_VERBOSE
	pargeo::timer t; t.start();
#endif

	pargeo::kdNode<3, pointT>* tree =
	  buildKdTree<3, pointT>(P, true, 128); // todo manual leaf size

#ifdef SAMPLE_HULL_VERBOSE
	std::cout << " build-tree-time = " << t.get_next() << "\n";
#endif

	parlay::sequence<std::atomic<bool>> flag(P.size());
	parlay::parallel_for(0, flag.size(), [&](size_t i){flag[i] = false;});

	auto facets = parlay::tabulate(F.size(), [&](size_t i) {
					   return facet(F[i].a, F[i].b, F[i].c);
					 });

	auto interiorPt = (facets[0].a + facets[facets.size()-1].a) / 2;

	auto data = &P[0];
	parlay::parallel_for(0, F.size(), [&](size_t i) {
					    flagVisible<pointT>(tree, &facets[i], make_slice(flag), data);
				  });

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
