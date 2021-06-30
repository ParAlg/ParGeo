#pragma once

#include <math.h>
#include "pargeo/kdTree.h"
#include "pargeo/point.h"
#include "convexHull3d/facet.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/hash_table.h"

namespace pargeo {
  namespace hullInternal {

    parlay::sequence<fpoint<3>>
    filter1(parlay::sequence<fpoint<3>> &P, parlay::sequence<facet3d<fpoint<3>>> &F) {
      using namespace parlay;

      auto interiorPt = (F[0].a + F[F.size()-1].a) / 2;

      sequence<fpoint<3>> area(F.size());
      parallel_for(0, F.size(),
		   [&](size_t i){
		     auto f = F[i];
		     area[i] = pargeo::crossProduct3d(f.b-f.a, f.c-f.a);
		   });

      auto isOut = [&](fpoint<3> p) {
		     for (size_t i = 0; i < F.size(); ++i) {
		       auto f = F[i];
		       if (((f.a - interiorPt) - (p - interiorPt)).dot(area[i]) > 0) {
			 return true;
		       }
		     }
		     return false;
		   };

      return parlay::filter(make_slice(P), isOut);
    }

    template<class pt>
    struct facet : public facet3d<pt> {
    public:

      using floatT = typename pt::floatT;
      using T = facet3d<pt>;

      pt area;

      floatT getVolume(pt p) {
	return (T::a - p).dot(area);
      }

      facet(pt _a, pt _b, pt _c): T(_a, _b, _c) {
	area = crossProduct3d(T::b - T::a, T::c - T::a);
      }

    };

    std::tuple<float, float> nodeVol(kdNode<3, fpoint<3>>* r, facet<fpoint<3>>* f) {
      using floatT = typename fpoint<3>::floatT;

      floatT vMin = std::numeric_limits<floatT>::max();
      floatT vMax = -1;
      floatT vol;
      fpoint<3> c;

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

    void flagVisible(kdNode<3, fpoint<3>>* r,
		     facet<fpoint<3>>* f,
		     parlay::slice<bool*, bool*> flag,
		     fpoint<3>* base) {
      using floatT = typename fpoint<3>::floatT;

      auto vols = nodeVol(r, f);
      floatT volMin = std::get<0>(vols);
      floatT volMax = std::get<1>(vols);

      if (volMin <= 0 && volMax <= 0) {
	return;
      } else if (volMin <= 0 && volMax > 0) {
	// not intersect
	if (!r->isLeaf()) {
	  flagVisible(r->L(), f, flag, base);
	  flagVisible(r->R(), f, flag, base);
	} else {
	  for (size_t i = 0; i < r->size(); ++ i) {
	    auto pVol = f->getVolume(*r->at(i));
	    if (pVol >= 0 ||
		f->a == *r->at(i) ||
		f->b == *r->at(i) ||
		f->c == *r->at(i)
		) {
	      flag[r->at(i) - base] = true;
	    }
	  }
	}
      } else {
	// node visible
	for (size_t i = 0; i < r->size(); ++ i) {
	  auto pVol = f->getVolume(*r->at(i));
	  if (pVol >= 0 ||
	      f->a == *r->at(i) ||
	      f->b == *r->at(i) ||
	      f->c == *r->at(i)
	      ) {
	    flag[r->at(i) - base] = true;
	  }
	}
      }

    }

    parlay::sequence<fpoint<3>>
    filter2(parlay::sequence<fpoint<3>> &P, parlay::sequence<facet3d<fpoint<3>>> &F) {
      using namespace parlay;

      kdNode<3, fpoint<3>>* tree =
	buildKdt<3, fpoint<3>>(P, true, false); // todo manual leaf size

      sequence<bool> flag = tabulate(P.size(), [&](size_t i) {
	  return false;
	});

      auto facets = tabulate(F.size(), [&](size_t i) {
	  return facet(F[i].a, F[i].b, F[i].c);
	});

      auto interiorPt = (facets[0].a + facets[facets.size()-1].a) / 2;

      parallel_for(0, F.size(), [&](size_t i) {
	  flagVisible(tree, &facets[i], make_slice(flag), P.data());
	});

      return parlay::pack(P, flag);
    }

  } // End namespace hullInternal
} // End namespace pargeo
