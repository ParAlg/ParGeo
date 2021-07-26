#include <atomic>

#include "pargeo/point.h"
#include "pargeo/getTime.h"
#include "pargeo/parlayAddon.h"

#include "parlay/sequence.h"

#include "convexHull2d/randInc/hull.h"

// #define RAND_HULL_2D_VERBOSE

using namespace std;
using namespace parlay;
using namespace pargeo;

template<class pointT> class facet;

template<class pointT>
class vertex: public pointT {
public:
  using floatT = typename pointT::floatT;

  facet<vertex>* seeFacet;

  vertex(floatT* p): pointT(p) {};
  vertex(floatT* p, typename pointT::attT a): pointT(p, a) {};
  vertex(): pointT() {};
};

template<class vertex>
class facet {

private:

  inline typename vertex::floatT triArea(vertex a, vertex b, vertex c) {
    auto cross = [&](vertex p1, vertex p2) {
      return p1[0] * p2[1] - p1[1] * p2[0];
    };
    return cross((b-a), (c-a));
  }

public:
  /* reservation */

  std::atomic<size_t> token;

  void reset() {
    token = -1;
  }

  template<class T>
  void reserve(T* pt) {
    parlay::write_min(&token, size_t(pt), std::less<size_t>());
  }

  template<class T>
  bool isReserved(T* pt) {
    return size_t(pt) == token;
  }

  /* facet pointers */
  facet* next; facet* prev;

  /* vertices (p1->p2) is clockwise */
  vertex* p1; vertex* p2;

  sequence<vertex*> seeList;

  facet(vertex* _p1, vertex* _p2): p1(_p1), p2(_p2) {
    reset();
  }

  bool visibleFrom(vertex p) {
    return triArea(*p1, *p2, p) > p.eps;
  }

};

template<class pointT>
pair<facet<vertex<pointT>>*, facet<vertex<pointT>>*>
findVisible(facet<vertex<pointT>>* head, vertex<pointT> p) {
  using facet = facet<vertex<pointT>>;

  while (head->prev->visibleFrom(p)) head = head->prev;

  auto ptr = head;
  facet* start = NULL;
  size_t n = 0;
  do {
    if (n <= 0) {
      if (ptr->visibleFrom(p)) {
        n++;
        start = ptr;
      }
    } else {
      if (!ptr->visibleFrom(p)) {
        return make_pair(start, ptr);
      }
      n++;
    }
    ptr = ptr->next;
  } while (ptr != head);
  return make_pair((facet*)NULL, (facet*)NULL);
}

template<class pointT>
parlay::sequence<size_t>
pargeo::hull2d::randInc::compute(parlay::slice<pointT*, pointT*> P, size_t numProc) {

  using vertex = vertex<pointT>;
  using facet = facet<vertex>;

#ifdef RAND_HULL_2D_VERBOSE
  timer t; t.start();
#endif

  /* Initialize vertex array */

  size_t n = P.size();

  sequence<vertex> V = parlay::tabulate(n, [&](size_t i){
    return vertex(P[i].coords());
  });

  /* Find extrema and initialize a hull */

  auto xLt = [&] (vertex a, vertex b) {
    return (a[0] < b[0]) || ((a[0] == b[0]) && (a[1] < b[1]));};

  auto yLt = [&] (vertex a, vertex b) {
    return (a[1] < b[1]) || ((a[1] == b[1]) && (a[0] < b[0]));};

  auto xExt = minmax_element(V, xLt);
  size_t xMin = xExt.first - V.data();
  size_t xMax = xExt.second - V.data();

  auto yExt = minmax_element(V, yLt);
  size_t yMin = yExt.first - V.data();
  size_t yMax = yExt.second - V.data();

  /* Create initial hull structure */

  auto f0 = new facet(&V[xMin], &V[yMax]);
  auto f1 = new facet(&V[yMax], &V[xMax]);
  auto f2 = new facet(&V[xMax], &V[yMin]);
  auto f3 = new facet(&V[yMin], &V[xMin]);
  auto H = f0;

  f0->prev = f3; f0->next = f1;
  f1->prev = f0; f1->next = f2;
  f2->prev = f1; f2->next = f3;
  f3->prev = f2; f3->next = f0;

  /* Split points based on the facets */

  // mark which facet the point belongs, 0, 1, 2, 3
  sequence<size_t> flag = parlay::tabulate(n, [&](size_t i) {
    if (f0->visibleFrom(V[i])) {
      V[i].seeFacet = f0;
      return (size_t)0;
    } else if (f1->visibleFrom(V[i])) {
      V[i].seeFacet = f1;
      return (size_t)1;
    } else if (f2->visibleFrom(V[i])) {
      V[i].seeFacet = f2;
      return (size_t)2;
    } else if (f3->visibleFrom(V[i])) {
      V[i].seeFacet = f3;
      return (size_t)3;
    } else { // in hull, sees no facet
      V[i].seeFacet = nullptr;
      return (size_t)4;
    }
  });

  sequence<vertex*> pointers =
    parlay::tabulate(n, [&](size_t i) {
      return &V[i];});

  sequence<sequence<vertex*>> chunks =
    parlay::split_k(4, pointers, flag);

#ifdef RAND_HULL_2D_VERBOSE
  std::cout << "first-slit = ";
  for (auto c: chunks) {
    std::cout << c.size() << ", ";
  }
  std::cout << "\n";
#endif

  f0->seeList = std::move(chunks[0]);
  f1->seeList = std::move(chunks[1]);
  f2->seeList = std::move(chunks[2]);
  f3->seeList = std::move(chunks[3]);

#ifdef RAND_HULL_2D_VERBOSE
  std::cout << "x-extrema = " << xMin << ", " << xMax << "\n";
  std::cout << "y-extrema = " << yMin << ", " << yMax << "\n";
  std::cout << "init-time = " << t.get_next() << "\n";
#endif

  // use a mapping to keep track of unprocessed nodes
  // and pack at the end of every round
  sequence<vertex*> mapping = std::move(pointers);
  size_t remaining = n;

  /* A packing function for pointers to the remaining points
   */
  auto packRemain = [&]() {

    auto flag = parlay::tabulate(remaining+1, [&](size_t i){
      return 0;});

    parallel_for(0, remaining, [&](size_t i) {
      if (mapping[i]->seeFacet != nullptr)
	flag[i] = 1;});

    size_t M2 = scan_inplace(make_slice(flag), parlay::addm<size_t>());
    flag[remaining] = M2;

    if (M2 <= 0) {
      remaining = M2;
      return;}

    auto map2 = parlay::sequence<vertex*>(M2);
    parallel_for(0, remaining,
		 [&](size_t i){
		   if (flag[i] != flag[i+1]) map2[flag[i]] = mapping[i];
		 });

    mapping = std::move(map2);
    remaining = M2;
  };

  packRemain();

  if (!numProc) numProc = parlay::num_workers();
  size_t batch = numProc * 4;

  while (remaining > 0) {

#ifdef RAND_HULL_2D_VERBOSE
    std::cout << "remaining = " << remaining << "\n";
#endif

    /* each interation processes a batch */

    size_t b = min(batch, remaining);

    /* reservation */

    parallel_for(0, b, [&](size_t i) {
      auto facets = findVisible(mapping[i]->seeFacet, *mapping[i]);
      facet* start = facets.first;
      facet*   end = facets.second;

      if (start && end) {
	auto ptr = start->prev;
	do {
	  ptr->reserve(mapping[i]);
	  ptr = ptr->next;
	} while (ptr != end->next);
      }
    });

    sequence<vertex*> succ(b);

    /* check reservation */

    parallel_for(0, b, [&](size_t i) {
      auto facets = findVisible(mapping[i]->seeFacet, *mapping[i]);
      facet* start = facets.first;
      facet*   end = facets.second;

      if (start && end) {
	auto ptr = start->prev;
	bool ok = true;
	do {
	  if (!ptr->isReserved(mapping[i])) {
	    ok = false;
	  } else {
	    ptr->reset();
	  }
	  ptr = ptr->next;
	} while (ptr != end->next);

	if (ok) succ[i] = mapping[i];
	else succ[i] = nullptr;
      }
    });

    succ = parlay::filter(make_slice(succ), [&](vertex* v){return v!=nullptr;});

#ifdef RAND_HULL_2D_VERBOSE
    std::cout << "succ = " << succ.size() << "/" << b << "\n";
#endif

    /* process successful */

    parallel_for(0, succ.size(), [&](size_t i) {
      auto facets = findVisible(succ[i]->seeFacet, *succ[i]);
      facet* start = facets.first;
      facet*   end = facets.second;

      if (start && end) {

	facet* new1 = new facet(start->p1, succ[i]);
	facet* new2 = new facet(succ[i], end->p1);

	// nullify seeFacet so succ[i] does not need to be processed again
	succ[i]->seeFacet = nullptr;

	// reservation = [start->prev, end]
	// the actual processing range is [start, end->prev]

	// gather necessary point arrays of the facets

	sequence<sequence<vertex*>> Q0;
	auto ptr = start;
	do {
	  Q0.push_back(std::move(ptr->seeList));
	  ptr = ptr->next;
	} while (ptr != end);
	sequence<vertex*> Q = parlay::flatten(make_slice(Q0));

	// redistribution

	sequence<size_t> flag =
	  parlay::tabulate(Q.size(), [&](size_t i) {
	    if (new1->visibleFrom(*Q[i])) {
	      Q[i]->seeFacet = new1;
	      return (size_t)0;
	    } else if (new2->visibleFrom(*Q[i])) {
	      Q[i]->seeFacet = new2;
	      return (size_t)1;
	    } else {
	      Q[i]->seeFacet = nullptr;
	      return (size_t)2;
	    }
	  });

	auto chunks = parlay::split_k(3, Q, flag);
	new1->seeList = std::move(chunks[0]);
	new2->seeList = std::move(chunks[1]);

	/* update hull*/

	start->prev->next = new1; new1->prev = start->prev;
	new1->next = new2; new2->prev = new1;
	new2->next = end; end->prev = new2;
	H = new1;

	ptr = start;
	do {
	  auto next = ptr->next;
	  delete ptr;
	  ptr = next;
	} while (ptr != end);
      }
    });

    packRemain();

  }

#ifdef RAND_HULL_2D_VERBOSE
  /* calculate hull size */

  size_t h = 0;
  auto ptr = H;
  do {
    h++;
    ptr = ptr->next;
  } while (ptr != H);

  std::cout << "hull-size = " << h << "\n";
#endif

  return parlay::sequence<size_t>(); // todo note this is dummy return
}

template parlay::sequence<size_t>
pargeo::hull2d::randInc::compute(parlay::slice<pargeo::fpoint<2>*, pargeo::fpoint<2>*>, size_t);

template parlay::sequence<size_t>
pargeo::hull2d::randInc::compute(parlay::slice<pargeo::point<2>*, pargeo::point<2>*>, size_t);
