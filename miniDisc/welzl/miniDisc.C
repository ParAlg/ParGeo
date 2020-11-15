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

#include <vector>
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "miniDisc.h"
#include "geometry.h"
#include "check.h"

using namespace std;

template<int dim>
ball<dim> miniDiscImpl(point<dim>* P, intT n, vector<point<dim>>& support, ball<dim> B) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  if (B.isEmpty()) {
    if (support.size() == 0) {
      B = ballT(P, 2);
    } else if (support.size() == 1) {
      support.push_back(P[0]);
      B = ballT(&support[0], support.size());
      support.pop_back();
    } else { //>=2
      B = ballT(&support[0], support.size());
    }
  }

  if (B.size() == dim+1) {
    return B;
  }

  for (intT i=0; i<n; ++i) {
    //cout << "i = " << i << endl;
    if (!B.contain(P[i])) {
      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();
      support.push_back(P[i]);
      B = miniDiscImpl(P, i, support, B);
      support.pop_back();
      //todo move to front
      //swap(P[dim-support.size()], P[i]);
      //swap(P[i], P[i-1]);
    }
  }

  return B;
}

template<int dim>
void miniDisc(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef circle discT;
  typedef ball<dim> ballT;

  //static const bool serial = true;
  static const bool noRandom = true;
  //static const bool moveToFront = true;

  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;

  timing t0;t0.start();
  if(!noRandom) {
    cout << "permuting points" << endl;
    randPerm(P, n);
  }

  auto support = vector<pointT>();
  auto D = miniDiscImpl(P, n, support, ballT());
  cout << D.radius() << ", center = " << D.center() << endl;

  cout << "total-time = " << t0.stop() << endl;
  cout << endl;

  check<dim,ballT>(&D, P, n);
}

template void miniDisc<2>(point<2>*, intT);
template void miniDisc<3>(point<3>*, intT);
template void miniDisc<4>(point<4>*, intT);
template void miniDisc<5>(point<5>*, intT);
template void miniDisc<6>(point<6>*, intT);
template void miniDisc<7>(point<7>*, intT);
template void miniDisc<8>(point<8>*, intT);
template void miniDisc<9>(point<9>*, intT);
