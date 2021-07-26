// This code is part of the Pargeo Library
// Copyright (c) 2021 Yiqiu Wang and the Pargeo Team
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

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "convexHull3d/parallelQuickHull/hullImpl.h"
#include "convexHull3d/parallelQuickHull/linkedFacet.h"

// #define HULL_PARALLEL_VERBOSE

namespace pargeo {
  namespace hull3d {
    namespace parallelQuickHull {

template<class linkedFacet3d, class vertex3d, class pointT>
void quickHullParallel(pargeo::hull3d::parallelQuickHull::hullTopology<pointT> *context,
		       size_t numProc = 0,
		       bool randomized = false) {
  using namespace pargeo;
  using namespace parlay;

  if (numProc == 0)
    numProc = parlay::num_workers();

#ifdef HULL_PARALLEL_VERBOSE
  pargeo::timer t;
  size_t round = 0;
  size_t serRound = 0;
  pargeo::timer tr;
  double serTime = 0;
  double parTime = 0;
  double apexTime = 0;
  double frontierTime = 0;
  double splitTime = 0;
  double totalApex = 0;
  double succApex = 0;
#endif

  while (true) {

#ifdef HULL_PARALLEL_VERBOSE
    //context->stats();
    round ++;
#endif

    // todo find better threshold
    if (context->hullSize() < 64) {

#ifdef HULL_PARALLEL_VERBOSE
      serRound ++;
      tr.start();
      t.start();
#endif

    loopStart:
      linkedFacet3d* f0 = context->hullSize() < 512 ? context->facetWalk() : nullptr;
      vertex3d *apex;
      if (randomized)
	apex = context->randomApex(f0);
      else
	apex = context->furthestApexParallel(f0);

#ifdef HULL_PARALLEL_VERBOSE
      apexTime += t.get_next();
#endif

      if (!apex) break;

      auto frontier = context->computeFrontier(*apex);
      auto frontierEdges = std::move(std::get<0>(frontier));
      auto facetsBeneath = std::move(std::get<1>(frontier));

#ifdef HULL_PARALLEL_VERBOSE
      frontierTime += t.get_next();
#endif

      // Check for frontier error, usually caused by numerical errors in rare cases
      for(size_t i=0; i<frontierEdges.size(); ++i) {
	auto nv = frontierEdges.at((i+1)%frontierEdges.size());
	auto cv = frontierEdges.at(i);
	if (cv.b != nv.a) {
	  apex->seeFacet->clear();
	  apex->seeFacet = nullptr;
	  goto loopStart;
	}
      }

      // Create new facets
      auto newFacets = sequence<linkedFacet3d*>(frontierEdges.size());

      for (size_t i=0; i<frontierEdges.size(); ++i) {
	typename pargeo::hull3d::parallelQuickHull::hullTopology<pointT>::_edge e = frontierEdges.at(i);
	newFacets[i] = new linkedFacet3d(e.a, e.b, *apex);
      }

      // Connect new facets
      for (size_t i=0; i<frontierEdges.size(); ++i) {
	context->linkFacet(newFacets[i],
		  newFacets[(i+1)%frontierEdges.size()],
		  frontierEdges.at(i).ff,
		  newFacets[(i-1+frontierEdges.size())%frontierEdges.size()]
		  );
      }

      context->setHull(newFacets[0]);

      context->redistributeParallel(make_slice(facetsBeneath), make_slice(newFacets));

#ifdef HULL_PARALLEL_VERBOSE
      splitTime += t.stop();
#endif

      // Delete existing facets
      for(int j=0; j<facetsBeneath.size(); ++j)
	delete facetsBeneath.at(j);

#ifdef HULL_PARALLEL_VERBOSE
      serTime += tr.stop();
#endif

    } else { // serial vs parallel if

#ifdef HULL_PARALLEL_VERBOSE
      tr.start();
      t.start();
#endif

      // the apex collecting function should also return how difficult it was to
      // collect the apexes, if it's very difficult, maybe should pack the vertices
      // and find them instead
      if (context->updateMap() <= 0) break;
      sequence<vertex3d*> apexes0;
      if (randomized)
	apexes0 = std::move(context->randomApexes(numProc * 8)); // todo tune
      else
	apexes0 = std::move(context->furthestApexesParallel3(numProc * 8)); // todo tune
      // sequence<vertex3d> apexes0 = context->furthestApexes(numProc * 8); // todo tune
      // sequence<vertex3d> apexes0 = context->randomApexes(numProc*8); // not as good
      //sequence<vertex3d> apexes0 = context->furthestApexesWithSkip(numProc*8); // not as good

#ifdef HULL_PARALLEL_VERBOSE
      apexTime += t.get_next();
#endif
      if (apexes0.size() <= 0) break;

      size_t numApex0 = apexes0.size();
      sequence<sequence<typename hull3d::parallelQuickHull::hullTopology<pointT>::_edge>> FE0(numApex0);
      sequence<sequence<linkedFacet3d*>> FB0(numApex0);

      // Compute frontier and reserve facets
      parallel_for(0, numApex0, [&](size_t a) {
				 auto frontier = context->computeFrontierAndReserve(*apexes0[a]);

				 sequence<typename hull3d::parallelQuickHull::hullTopology<pointT>::_edge> frontierEdges = std::move(std::get<0>(frontier));
				 sequence<linkedFacet3d*> facetsBeneath = std::move(std::get<1>(frontier));
				 FE0[a] = std::move(frontierEdges);
				 FB0[a] = std::move(facetsBeneath);
			       });
#ifdef HULL_PARALLEL_VERBOSE
      frontierTime += t.get_next();
#endif

      // Check reservation for success, then reset reservation
      // Also check for numerical errors, set to unsuccessful if detected
      sequence<bool> success(numApex0);
      parallel_for(0, numApex0, [&](size_t a) {

				 // Check reservation
				 if (!context->confirmReservation( *apexes0[a], FB0[a].cut(0, FB0[a].size()) )) {
				   success[a] = false;
				 } else {
				   success[a] = true;
				 }

				 // Check for numerical errors in the frontier
				 for(size_t i=0; i<FE0[a].size(); ++i) {
				   auto nv = FE0[a].at((i+1)%FE0[a].size());
				   auto cv = FE0[a].at(i);
				   if (cv.b != nv.a) {
				     // Error found in point a
				     linkedFacet3d* f = apexes0[a]->seeFacet;
				     parallel_for(0, f->numPts(),
						  [&](size_t i){
						    f->pts(i)->seeFacet = nullptr;
						  });
				     f->clear();
				     apexes0[a]->seeFacet = nullptr;
				     success[a] = false;
				     break;
				   }
				 }

				});
      parallel_for(0, numApex0, [&](size_t a) {
				 context->resetReservation(*apexes0[a], FB0[a].cut(0, FB0[a].size()));
			       });

      auto apexes = parlay::pack(make_slice(apexes0), make_slice(success));
      auto FE = parlay::pack(make_slice(FE0), make_slice(success));
      auto FB = parlay::pack(make_slice(FB0), make_slice(success));
      size_t numApex = apexes.size();

#ifdef HULL_PARALLEL_VERBOSE
      //std::cout << "reservation = " << numApex << " / " << apexes0.size() << "\n";
      totalApex += apexes0.size();
      succApex += numApex;
#endif

      sequence<int> increase(numApex, 0);

      // Process the successful points
      parallel_for(0, numApex, [&](size_t a) {
				   // Create new facets
				   auto newFacets = sequence<linkedFacet3d*>(FE[a].size());

				   for (size_t i=0; i<FE[a].size(); ++i) {
				     typename hull3d::parallelQuickHull::hullTopology<pointT>::_edge e = FE[a].at(i);
				     newFacets[i] = new linkedFacet3d(e.a, e.b, *apexes[a]);
				   }

				   // Connect new facets
				   for (size_t i=0; i<FE[a].size(); ++i) {
				     context->linkFacet(newFacets[i],
					       newFacets[(i+1)%FE[a].size()],
					       FE[a].at(i).ff,
					       newFacets[(i-1+FE[a].size())%FE[a].size()]
					       );
				   }

				   context->setHull(newFacets[0]); // todo data race

				   context->redistributeParallel(FB[a].cut(0, FB[a].size()), make_slice(newFacets));

				   increase[a] = FE[a].size() - FB[a].size();
			       });

      // if (!context->checkReset()) {
      // 	throw std::runtime_error("Facet not reset properly");
      // }

#ifdef HULL_PARALLEL_VERBOSE
      splitTime += t.stop();
#endif

      // Clean up

      parallel_for(0, numApex0, [&](size_t a) {
				 if (success[a]) {
				   // Delete existing facets
				   for(int j=0; j<FB0[a].size(); ++j) {
				     auto f = FB0[a].at(j);
				     delete FB0[a].at(j);
				   }
				 }
			       });

      // todo may not be necessary
      if (context->numConflict() <= 0) break;

#ifdef HULL_PARALLEL_VERBOSE
      parTime += tr.stop();
#endif

    } // End serial/parallel if

  }
#ifdef HULL_PARALLEL_VERBOSE
  std::cout << "------\n";
  std::cout << "apex-time = " << apexTime << "\n";
  std::cout << "frontier-time = " << frontierTime << "\n";
  std::cout << "create-split-time = " << splitTime << "\n";
  std::cout << "------\n";
  std::cout << "#-rounds = " << round << "\n";
  std::cout << "#-ser-rounds = " << serRound << "\n";
  std::cout << "ser-time = " << serTime << "\n";
  std::cout << "#-par-rounds = " << round - serRound << "\n";
  std::cout << "par-time = " << parTime << "\n";
  std::cout << "------\n";
  std::cout << "reservation-count = " << totalApex << "\n";
  std::cout << "reservation-succ-rate = " << succApex / totalApex << "\n";
  std::cout << "------\n";
  std::cout << "hull-size = " << context->hullSize() << "\n";
#endif

  //context->writeHull();
}

    }
  }
}
