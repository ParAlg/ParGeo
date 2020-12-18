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
#include "pbbs/gettime.h"
#include "pbbs/sequence.h"
#include "geometry.h"
#include "hull.h"
using namespace std;

struct pointNode;

struct facet {
  facet* next;
  facet* prev;
  point2d p1;
  point2d p2;//(p1->p2) is clockwise
  vector<pointNode*> seeList;//points that can see facet
  facet(point2d p11, point2d p22): p1(p11), p2(p22) {}
  bool visibleFrom(point2d p) {return triArea(p1, p2, p) > 1e-10;}//todo numerical stability needs better hack
};

static std::ostream& operator<<(std::ostream& os, const facet f) {
  os << "(" << f.p1 << ", " << f.p2 << ")";
  return os;
}

struct pointNode {
  point2d p;
  facet* seeFacet;//maintain one edge visible, change from time to time
  //don't have to maintain all, can search pretty easily
  pointNode(point2d pp): p(pp), seeFacet(NULL) {};
  pointNode() {};
};

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
  static bool brute = false;//visibility check
  static bool verify = false;
  static bool localPivot = false;

  timing t; t.start();

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
  if(verbose) {
    cout << "initial-hull = ";
    printHull(H, H);
  }

  pointNode* PN = newA(pointNode, n-3);
  for(intT i=0; i<n-3; ++i) {
    PN[i] = pointNode(P[i+3]);
    auto ptr = H;
    do {
      if (ptr->visibleFrom(P[i+3])) {
        PN[i].seeFacet = ptr;
        ptr->seeList.push_back(&PN[i]);//todo this vector needs to be simplified (mem management)
        break; //each point only records one visible facet
      }
      ptr = ptr->next;
    } while (ptr != H);
  }

  //given facet H, find furthest visible point (only among recorded)
  auto findPivot = [&](facet* H)
		   {
		     auto f = H;
		     while (f->seeList.size()<=0 && f->next!=H) f = f->next;
		     if (f->seeList.size() <= 0) return (pointNode*)NULL;
		     point2d l = f->p1;
		     point2d r = f->p2;
		     auto triangArea = [&](intT idx)
				       {
					 return triArea(l, r, f->seeList[idx]->p);
				       };
		     intT idx = sequence::maxIndex<double>((intT)0, (intT)f->seeList.size(), greater<floatT>(), triangArea);
		     return f->seeList[idx];
		   };

  intT i = 0;
  while(1) {
    if(verbose) cout << "--- iter" << i << endl;

    //find a pivot point to be processed next
    pointNode pr;
    if (localPivot) {
      pointNode* prp = findPivot(H);
      if (!prp) break;
      pr = *prp;
    } else {
      if (i>=n-3) break;
      pr = PN[i];
      if (!pr.seeFacet) {
	i++;
	continue;}
    }
    i++;

    if(verbose) cout << " pr = " << pr.p << endl;

    //find range of visible facets [left, right)
    pair<facet*, facet*> conflicts;
    if (brute) {
      conflicts = findVisible(H, pr.p);
    } else {
      conflicts = findVisible(pr.seeFacet, pr.p);
    }

    facet* start = conflicts.first;
    facet*   end = conflicts.second;
    if (start && end) {
      //compute new faces
      facet* new1 = new facet(start->p1, pr.p);
      facet* new2 = new facet(pr.p, end->p1);

      //go through list of facets pending to be deleted
      if (verbose) {
        cout << "deleting = ";
        auto ptr = start;
        do {
          cout << *ptr << ", ";
          ptr = ptr->next;
        } while (ptr != end);
        cout << endl;
      }

      //if not finding visible facets by bruteforce
      // update visibility pointers
      if (!brute) {
        auto ptr = start;
        do {
          //update points that see them
          for(intT j=0; j<ptr->seeList.size(); ++j) {
            auto seePt = ptr->seeList[j];
            if (seePt->seeFacet == ptr) seePt->seeFacet = NULL;

            if (new1->visibleFrom(seePt->p)) {
              if (!seePt->seeFacet) {
                seePt->seeFacet = new1;
              }
              new1->seeList.push_back(seePt);
            } else if (new2->visibleFrom(seePt->p)) {
              if (!seePt->seeFacet) {
                seePt->seeFacet = new2;
              }
              new2->seeList.push_back(seePt);
            }
          }

          ptr = ptr->next;
        } while (ptr != end);
      }

      if(verbose) cout << "adding = " << *new1 << ", " << *new2 << endl;

      //update hull
      start->prev->next = new1; new1->prev = start->prev;
      new1->next = new2; new2->prev = new1;
      new2->next = end; end->prev = new2;
      H = new1;
      delHull(start, end);

      if(verbose) {
        cout << "hull = ";
        printHull(H, H);
      }
    }
  }

  free(PN);
  cout << "hull-time = " << t.next() << endl;

  if (verify) {
    for (intT i=0; i<n; ++i) {
      auto ptr = H;
      do {
        if (ptr->visibleFrom(P[i])) {
          cout << "wrong hull, " << P[i] << " visible from " << *ptr << endl;
          cout << "triArea = " << triArea(ptr->p1, ptr->p2, P[i]) << endl;
          abort();
        }
        ptr = ptr->next;
      } while (ptr != H);
    }
    cout << "hull verified, time = " << t.stop() << endl;
  }

  intT hSize = 0;
  auto ptr = H;
  do {
    hSize ++;
    ptr = ptr->next;
  } while (ptr != H);
  cout << "hull size = " << hSize << endl;

  intT* I = newA(intT, n);
  return _seq<intT>(I,1);//dummy
}
