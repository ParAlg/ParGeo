// This code is part of the Pargeo Library
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

//#define WRITE // Write to file, visualize using python3 plot.py
//#define VERBOSE
#define SILENT

#ifdef WRITE
#include <iostream>
#include <fstream>
#endif

#include "parlay/sequence.h"
#include "geometry/point.h"
#include "common/get_time.h"
#include "hull.h"
#include "split.h"
#include "hullTopology.h"

using namespace parlay;
using namespace parlay::internal;

template<class linkedFacet3d, class vertex3d, class origin3d>
void incrementHull3dSerial(_hull<linkedFacet3d, vertex3d, origin3d> *context) {

  timer t;// t.start();

  size_t errors = 0;
  size_t round = 0;
  double apexTime = 0;
  double frontierTime = 0;
  double createTime = 0;
  double splitTime = 0;

  while (true) {
    t.start();

  loopStart:
    // Serial

#ifdef WRITE
    context->writeHull();
#endif

    //vertex3d apex = context->randomApex();
    vertex3d apex = context->furthestApex();

#ifdef VERBOSE
    cout << ">>>>>>>>>" << endl;
    //context->printHull();
    cout << "apex chosen = " << apex.attribute.i << endl;
#endif
    apexTime += t.get_next();

    if (apex.isEmpty()) break;

    round ++;

    auto frontier = context->computeFrontier(apex);
    auto frontierEdges = get<0>(frontier);
    auto facetsBeneath = get<1>(frontier);

    frontierTime += t.get_next();

    for(size_t i=0; i<frontierEdges->size(); ++i) {
      auto nv = frontierEdges->at((i+1)%frontierEdges->size());
      auto cv = frontierEdges->at(i);
      if (cv.b != nv.a) {
	apex.attribute.seeFacet->clear();
	errors ++;
	goto loopStart;
      }
    }

#ifdef VERBOSE
    cout << "visible = " << *apex.attribute.seeFacet << ": ";
    for (auto x: *apex.attribute.seeFacet->seeList) {
      cout << x.attribute.i << " ";
    }
    cout << endl;

    cout << "frontier = ";
    for(auto e: *frontierEdges) {
      cout << e.a.attribute.i << "," << e.b.attribute.i << " ";
    }
    cout << endl;

    cout << "to delete = ";
    for(auto f: *facetsBeneath)
      cout << *f << " ";
    cout << endl;
#endif

    // Create new facets
    auto newFacets = sequence<linkedFacet3d*>(frontierEdges->size());

    for (size_t i=0; i<frontierEdges->size(); ++i) {
      typename _hull<linkedFacet3d, vertex3d, origin3d>::_edge e = frontierEdges->at(i);
      newFacets[i] = new linkedFacet3d(e.a, e.b, apex);
    }

#ifdef VERBOSE
    cout << "to create = ";
    for (size_t i=0; i<frontierEdges->size(); ++i)
      cout << *(newFacets[i]) << " ";
    cout << endl;
#endif

    // Connect new facets
    for (size_t i=0; i<frontierEdges->size(); ++i) {
      context->linkFacet(newFacets[i],
		newFacets[(i+1)%frontierEdges->size()],
		frontierEdges->at(i).ff,
		newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		);
    }

    context->setStart(newFacets[0]);

    createTime += t.get_next();

    context->redistribute(facetsBeneath->cut(0, facetsBeneath->size()), make_slice(newFacets));

    splitTime += t.stop();

    // Delete existing facets
    for(int j=0; j<facetsBeneath->size(); ++j)
      delete facetsBeneath->at(j);

    delete frontierEdges;
    delete facetsBeneath;

  }

#ifndef SILENT
  cout << "init-time = " << initTime << endl;
  cout << "apex-time = " << apexTime << endl;
  cout << "frontier-time = " << frontierTime << endl;
  cout << "create-time = " << createTime << endl;
  cout << "split-time = " << splitTime << endl;
  cout << "#-rounds = " << round << endl;
  cout << "#-errors = " << errors << endl;
#endif

#ifdef VERBOSE
  cout << "hull-size = " << context->hullSizeDfs() << endl;
#endif

#ifdef WRITE
  context->writeHull();
#endif

}

#ifndef SERIAL

template<class linkedFacet3d, class vertex3d, class origin3d>
void incrementHull3d(_hull<linkedFacet3d, vertex3d, origin3d> *context, size_t numProc = 0) {

  if (numProc == 0)
    numProc = parlay::num_workers();

  timer t; t.start();

  double initTime = t.stop();

  size_t round = 0;
  size_t serRound = 0;
  timer tr;
  double serTime = 0;
  double parTime = 0;

  double apexTime = 0;
  double frontierTime = 0;
  double splitTime = 0;

  while (true) {
    round ++;

    if (context->hullSize() < 64) {

      serRound ++;
      tr.start();
      t.start();

    loopStart:
      linkedFacet3d* f0 = context->hullSize() < 512 ? context->facetWalk() : nullptr;
      vertex3d apex = context->furthestApex(f0);

      apexTime += t.get_next();

      if (apex.isEmpty()) break;

      auto frontier = context->computeFrontier(apex);
      auto frontierEdges = get<0>(frontier);
      auto facetsBeneath = get<1>(frontier);

      frontierTime += t.get_next();

      // Check for frontier error, usually caused by numerical errors in rare cases
      for(size_t i=0; i<frontierEdges->size(); ++i) {
	auto nv = frontierEdges->at((i+1)%frontierEdges->size());
	auto cv = frontierEdges->at(i);
	if (cv.b != nv.a) {
	  apex.attribute.seeFacet->clear();
	  goto loopStart;
	}
      }

      // Create new facets
      auto newFacets = sequence<linkedFacet3d*>(frontierEdges->size());

      for (size_t i=0; i<frontierEdges->size(); ++i) {
	typename _hull<linkedFacet3d, vertex3d, origin3d>::_edge e = frontierEdges->at(i);
	newFacets[i] = new linkedFacet3d(e.a, e.b, apex);
      }

      // Connect new facets
      for (size_t i=0; i<frontierEdges->size(); ++i) {
	context->linkFacet(newFacets[i],
		  newFacets[(i+1)%frontierEdges->size()],
		  frontierEdges->at(i).ff,
		  newFacets[(i-1+frontierEdges->size())%frontierEdges->size()]
		  );
      }

      context->setStart(newFacets[0]);

      context->redistribute(facetsBeneath->cut(0, facetsBeneath->size()), make_slice(newFacets));

      splitTime += t.stop();

      // Delete existing facets
      for(int j=0; j<facetsBeneath->size(); ++j)
	delete facetsBeneath->at(j);

      delete frontierEdges;
      delete facetsBeneath;

      serTime += tr.stop();
    } else {
      tr.start();
      t.start();

      // todo points chosen are close, which is bad
      sequence<vertex3d> apexes0 = context->furthestApexes(numProc); // todo tune
      //sequence<vertex3d> apexes0 = cg->randomApexes(context->hullSize()); // todo tune
#ifdef VERBOSE
      size_t np = apexes0.size();
#endif

      apexTime += t.get_next();

      if (apexes0.size() <= 0) break;

      size_t numApex0 = apexes0.size();
      sequence<sequence<typename _hull<linkedFacet3d, vertex3d, origin3d>::_edge>*> FE0(numApex0);
      sequence<sequence<linkedFacet3d*>*> FB0(numApex0);

      // Compute frontier and reserve facets
      parallel_for(0, numApex0, [&](size_t a) {
				 auto frontier = context->computeFrontierAndReserve(apexes0[a]);

				 sequence<typename _hull<linkedFacet3d, vertex3d, origin3d>::_edge>* frontierEdges = get<0>(frontier);
				 sequence<linkedFacet3d*>* facetsBeneath = get<1>(frontier);
				 FE0[a] = frontierEdges;
				 FB0[a] = facetsBeneath;
			       });

      frontierTime += t.get_next();

      // Check reservation for success, then reset reservation
      // Also check for numerical errors, set to unsuccessful if detected
      sequence<bool> success(numApex0);
      parallel_for(0, numApex0, [&](size_t a) {

				 // Check reservation
				 if (!context->confirmReservation( apexes0[a], FB0[a]->cut(0, FB0[a]->size()) )) {
				   success[a] = false;
				 } else {
				   success[a] = true;
				 }

				 // Check for numerical errors in the frontier
				 for(size_t i=0; i<FE0[a]->size(); ++i) {
				   auto nv = FE0[a]->at((i+1)%FE0[a]->size());
				   auto cv = FE0[a]->at(i);
				   if (cv.b != nv.a) {
				     apexes0[a].attribute.seeFacet->clear();
				     success[a] = false;
				   }
				 }

				});
      parallel_for(0, numApex0, [&](size_t a) {
				 context->resetReservation(apexes0[a], FB0[a]->cut(0, FB0[a]->size()));
			       });

      auto apexes = parlay::pack(make_slice(apexes0), success);
      auto FE = parlay::pack(make_slice(FE0), success);
      auto FB = parlay::pack(make_slice(FB0), success);
      size_t numApex = apexes.size();
#ifdef VERBOSE
      cout << numApex << "/" << np << endl;
#endif
      sequence<int> increase(numApex, 0);

      // Process the successful points
      parallel_for(0, numApex, [&](size_t a) {
				   // Create new facets
				   auto newFacets = sequence<linkedFacet3d*>(FE[a]->size());

				   for (size_t i=0; i<FE[a]->size(); ++i) {
				     typename _hull<linkedFacet3d, vertex3d, origin3d>::_edge e = FE[a]->at(i);
				     newFacets[i] = new linkedFacet3d(e.a, e.b, apexes[a]);
				   }

				   // Connect new facets
				   for (size_t i=0; i<FE[a]->size(); ++i) {
				     context->linkFacet(newFacets[i],
					       newFacets[(i+1)%FE[a]->size()],
					       FE[a]->at(i).ff,
					       newFacets[(i-1+FE[a]->size())%FE[a]->size()]
					       );
				   }

				   context->setStart(newFacets[0]); // todo data race

				   context->redistribute(FB[a]->cut(0, FB[a]->size()), make_slice(newFacets));

				   increase[a] = FE[a]->size() - FB[a]->size();
			       });

      //todo remove
      // if (!context->checkReset()) {// checking this is expensive O(nh)
      //   cout << "not reset " << endl; abort();
      // }

      splitTime += t.stop();

      // Clean up
      parallel_for(0, numApex0, [&](size_t a) {
				 if (success[a]) {
				   // Delete existing facets
				   for(int j=0; j<FB0[a]->size(); ++j)
				     delete FB0[a]->at(j);
				 }

				 delete FE0[a];
				 delete FB0[a];
			       });

      parTime += tr.stop();
    } // End serial/parallel if

  }
#ifndef SILENT
  cout << "init-time = " << initTime << endl;
  cout << "apex-time = " << apexTime << endl;
  cout << "frontier-time = " << frontierTime << endl;
  cout << "create-split-time = " << splitTime << endl;
  cout << "#-rounds = " << round << endl;
  cout << "#-ser-rounds = " << serRound << endl;
  cout << "ser-time = " << serTime << endl;
  cout << "#-par-rounds = " << round - serRound << endl;
  cout << "par-time = " << parTime << endl;
  cout << "hull-size = " << context->hullSize() << endl;
#endif

#ifdef WRITE
  context->writeHull();
#endif
}

#endif
