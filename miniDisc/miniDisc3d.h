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

#ifndef MINI_DISC_3D
#define MINI_DISC_3D

#include "pbbs/utils.h"
#include "pbbs/sequence.h"
#include "geometry.h"

template<int dim>
sphere<dim> miniDisc3DParallel(point<dim>* P, intT n) {
  return miniDisc3DSerial(P, n);//todo
}

template<int dim>
sphere<dim> miniDisc3DSerial(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef sphere<dim> sphereT;

  auto sphere = sphereT(P[0], P[1], P[2], P[3]);
  intT i = 4;

  while (i < n) {
    bool conflict = false;
    for (; i<n; ++i) {
      if (!sphere.contain(P[i])) {
        conflict = true;
        break;}
    }

    if (conflict) {
      cout << "conflict = " << i << endl;

      sphere = sphereT(P[0], P[1], P[2], P[i]);

      for (intT j=3; j<i; ++j) {//update 1
        if (!sphere.contain(P[j])) {
          sphere = sphereT(P[0], P[1], P[i], P[j]);

          for (intT k=2; k<j; ++k) {//update 2
            if (!sphere.contain(P[k])) {
              sphere = sphereT(P[0], P[i], P[j], P[k]);

              for (intT l=1; l<k; ++l) {//update 3
                if (!sphere.contain(P[l])) {
                  sphere = sphereT(P[i], P[j], P[k], P[l]);
                }}

            }}

        }}

    }
  }
  return sphere;
}

#endif
