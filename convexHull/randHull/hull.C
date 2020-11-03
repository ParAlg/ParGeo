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
#include <algorithm>
#include <iostream>
#include "pbbs/utils.h"
#include "geometry.h"
#include "hull.h"
using namespace std;

struct facet {
  facet* next;
  facet* prev;
  point2d p1;
  point2d p2;//(p1->p2) is clockwise
  facet(point2d p11, point2d p22): p1(p11), p2(p22) {}
  bool visibleFrom(point2d p) {return triArea(p1, p2, p) > 0.0;}
};

static std::ostream& operator<<(std::ostream& os, const facet f) {
  os << "(" << f.p1 << ", " << f.p2 << ")";
  return os;
}

pair<facet*, facet*> findVisible(facet* head, point2d p) {
  while (head->prev->visibleFrom(p)) head = head->prev;

  auto ptr = head;
  facet* start = NULL;
  intT n = 0;
  do {
    if (n <= 0) {
      if (ptr->visibleFrom(p)) {
        n++;
        start = ptr;
      }
    } else {
      if (!ptr->visibleFrom(p)) {
        //cout << "num conflicts = " << n << endl;
        return make_pair(start, ptr);
      }
      n++;
    }
    ptr = ptr->next;
  } while (ptr != head);
  //cout << "num conflicts = " << 0 << endl;
  return make_pair((facet*)NULL, (facet*)NULL);
}

void delHull(facet* start, facet* end) {
  auto ptr = start;
  do {
    auto tmp = ptr->next;
    delete ptr;
    ptr = tmp;
  } while (ptr != end);
}

void printHull(facet* start, facet* end) {
  auto ptr = start;
  do {
    cout << *ptr << " ";
    ptr = ptr->next;
  } while (ptr != end);
  cout << endl;
}

_seq<intT> hull(point2d* P, intT n) {
  static bool verbose = false;

  auto p0 = P[0].x() < P[1].x() ? P[0] : P[1];//left of p1
  auto p1 = P[0].x() < P[1].x() ? P[1] : P[0];
  auto p2 = P[2];

  facet* H;
  if (triArea(p0, p1, p2) > 0.0) {
    auto f0 = new facet(p0, p2);
    auto f1 = new facet(p2, p1);
    auto f2 = new facet(p1, p0);
    f0->next = f1; f1->prev = f0;
    f1->next = f2; f2->prev = f1;
    f2->next = f0; f0->prev = f2;
    H = f0;
  } else {
    auto f0 = new facet(p0, p1);
    auto f1 = new facet(p1, p2);
    auto f2 = new facet(p2, p0);
    f0->next = f1; f1->prev = f0;
    f1->next = f2; f2->prev = f1;
    f2->next = f0; f0->prev = f2;
    H = f0;
  }

  for (intT i=3; i<n; ++i) {
    if(verbose) cout << "---" << endl;

    auto pr = P[i];
    if(verbose) cout << i << " pr = " << pr << endl;

    auto conflicts = findVisible(H, pr);//bruteforce for now
    facet* start = conflicts.first;
    facet*   end = conflicts.second;

    if (start && end) {//has conflicts
      facet* new1 = new facet(start->p1, pr);
      facet* new2 = new facet(pr, end->p1);
      start->prev->next = new1; new1->prev = start->prev;
      new1->next = new2; new2->prev = new1;
      new2->next = end; end->prev = new2;
      H = new1;
      delHull(start, end);
      if(verbose) printHull(H, H);
    }
  }
  intT* I = newA(intT, n);
  return _seq<intT>(I,1);//dummy
}
