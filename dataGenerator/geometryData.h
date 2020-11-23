// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
//  (adapted from the Problem Based Benchmark Suite (PBBS)
//  by Guy Blelloch and the PBBS team)
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

#ifndef GEOMETRY_DATA_H
#define GEOMETRY_DATA_H

#include "geometry.h"
#include "dataGen.h"
#include "pbbs/parallel.h"

template<int dim>
point<dim> randNd(intT i, floatT scale=1) {
  uintT s[dim];
  s[0] = i;
  for (int j=1; j<dim; ++j) {
    s[j] = j*i + dataGen::hash<uintT>(s[j-1]);
  }
  floatT ss[dim];
  for (int j=0; j<dim; ++j) {
    ss[j] = 2*dataGen::hash<double>(s[j])-1;
  }
  return point<dim>(ss)*scale;
}

template<int dim>
point<dim> randInUnitSphere(intT i, floatT scale=1) {
  auto origin = point<dim>();
  for(int j=0; j<dim; ++j) origin[j] = 0;
  intT j = 0;
  point<dim> p;
  do {
    intT o = dataGen::hash<intT>(j++);
    p = randNd<dim>(o+i);
  } while (p.dist(origin) > 1.0);
  return p*scale;
}

template<int dim>
point<dim> randOnUnitSphere(intT i, floatT scale=1) {
  auto origin = point<dim>();
  for(int j=0; j<dim; ++j) origin[j] = 0;
  point<dim> v = randInUnitSphere<dim>(i);
  return (v/v.dist(origin))*scale;
}

point<2> randOnParabola(intT i, floatT scale=1) {
  point<2> p = randNd<2>(i);
  p[1] = p[0]*p[0];
  return p*scale;
}

point<3> randOnParabloid(intT i, floatT scale=1) {
  point<3> p = randNd<3>(i);
  p[2] = p[0]*p[0]+p[1]*p[1];
  return p*scale;
}

#endif
