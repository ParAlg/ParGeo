#include "parlay/sequence.h"
#include "pargeo/point.h"

#include "convexHull2d/bruteforce/hull.h"

template <class pointT>
inline typename pointT::floatT triArea(pointT a, pointT b, pointT c) {
  auto cross = [&](pointT p1, pointT p2) {
		 return p1[0] * p2[1] - p1[1] * p2[0];
	       };
  return cross((b-a), (c-a));
}

template<class pointT>
parlay::sequence<size_t>
pargeo::hull2d::bruteforce::compute(parlay::slice<pointT*, pointT*> P) {
  struct facet {
    size_t a;
    size_t b;

    facet(size_t _a, size_t _b):
      a(_a), b(_b) {}

    bool operator==(facet f2) {
      return (a == f2.a && b == f2.b) ||
	(b == f2.a && a == f2.b);
    }

    bool connect(facet f2) {
      return !(*this == f2) && (b == f2.a || b == f2.b);
    }

    facet getConnect(facet f2) {
      if (connect(f2)) {
	if (b == f2.a) return f2;
	else {
	  std::swap(f2.a, f2.b);
	  return f2;
	}
      } else {
	throw std::runtime_error("can't get disconnectet facet");
      }
    }
  };

  parlay::sequence<facet> H;

  for (size_t i = 0; i < P.size(); ++ i) {
    for (size_t j = i + 1; j < P.size(); ++ j) {
      bool dir;

      size_t l = 0;
      for (; l < P.size(); ++ l) {
	if (l == i || l == j) continue;
	dir = triArea(P[i], P[j], P[l]) > 0.0;
	break;
      }

      bool onHull = true;
      for (; l < P.size(); ++ l) {
	if (l == i || l == j) continue;
	if (triArea(P[i], P[j], P[l]) > 0.0 != dir) {
	  onHull = false;
	  break;
	}
      }

      if (onHull) H.emplace_back(i, j);
    }
  }

  parlay::sequence<size_t> H2;

  H2.push_back(H[0].a);
  facet prev = H[0];

  while (H2.size() < H.size()) {

    bool connected = false;

    for (auto f: H) {
      if (prev.connect(f)) {
	facet next = prev.getConnect(f);
	H2.push_back(next.a);
	prev = next;
	connected = true;
	break;
      }
    }

    if (!connected)
      throw std::runtime_error("hull is not connected");
  }

  return H2;
}

template parlay::sequence<size_t>
pargeo::hull2d::bruteforce::compute(parlay::slice<pargeo::fpoint<2>*, pargeo::fpoint<2>*>);

template parlay::sequence<size_t>
pargeo::hull2d::bruteforce::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*>);
