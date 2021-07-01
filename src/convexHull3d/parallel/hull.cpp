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

#include "convexHull3d/parallelHull.h"

#include "parlay/sequence.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"

#include "parallelHull.h"
#include "incremental.h"
#include "vertex.h"

parlay::sequence<pargeo::facet3d<pargeo::fpoint<3>>>
pargeo::hull3dParallel(parlay::sequence<pargeo::fpoint<3>> &P, size_t numProc) {
  using namespace std;
  using namespace parlay;
  using floatT = pargeo::fpoint<3>::floatT;
  using pointT = pargeo::fpoint<3>;
  using facetT = facet3d<pargeo::fpoint<3>>;
  using vertexT = pargeo::hullInternal::vertex;

  size_t n = P.size();

  sequence<vertexT> Q(P.size());
  parallel_for(0, P.size(), [&](size_t i) {
			      Q[i] = vertexT(P[i].coords());
			    });

  // Create an initial simplex
  auto origin = pointOrigin();

  auto linkedHull = new parallelHull<linkedFacet3d<vertexT>, vertexT, pointOrigin>(make_slice(Q), origin);

  incrementHull3d<linkedFacet3d<vertexT>, vertexT, pointOrigin>(linkedHull, numProc);

  // getHull will undo the translation of linkedHull
  auto out = sequence<facetT>();
  linkedHull->getHull<pointT>(out);

  delete linkedHull;
  return out;
}
