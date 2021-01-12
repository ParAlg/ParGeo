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
  // Data fields
  facet* next;
  facet* prev;
  point2d p1;
  point2d p2;//(p1->p2) is clockwise
  vector<pointNode*>* seeList;//points that can see facet

  // Constructors
  facet(point2d p11, point2d p22): p1(p11), p2(p22) {
    seeList = new vector<pointNode*>();}//todo allocation

  // Methods
  bool visibleFrom(point2d p) {return triArea(p1, p2, p) > numericKnob;}
  void push_back(pointNode* p) {seeList->push_back(p);}
  intT size() {return seeList->size();};
  pointNode* at(intT i) {return seeList->at(i);}
};

static std::ostream& operator<<(std::ostream& os, const facet f) {
  os << "(" << f.p1 << ", " << f.p2 << ")";
  return os;
}

struct pointNode {
  point2d p;
  facet* seeFacet;//maintain one edge visible, change from time to time
  //don't have to maintain all facets, can search pretty easily
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
        return make_pair(start, ptr);
      }
      n++;
    }
    ptr = ptr->next;
  } while (ptr != head);
  return make_pair((facet*)NULL, (facet*)NULL);
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
  static bool farPivot = true;

  auto facets = newA(facet, 2*n);

  auto newFacetLeft = [&](intT i, point2d p11, point2d p22) {
			facets[2*i] = facet(p11, p22);
			return &facets[2*i];
		      };
  auto newFacetRight = [&](intT i, point2d p11, point2d p22) {
			 facets[2*i+1] = facet(p11, p22);
			 return &facets[2*i+1];
		       };

  timing t; t.start();

  // Akl–Toussaint heuristic
  intT iTop, iBot, iLeft, iRight;
  floatT vTop, vBot, vLeft, vRight;
  vTop = floatMin();
  vBot = floatMax();
  vLeft = floatMax();
  vRight = floatMin();
  for(intT i=0; i<n; ++i) {
    floatT x = P[i].x();
    floatT y = P[i].y();
    if (y > vTop) {
      vTop = y;
      iTop = i;}
    if (y < vBot) {
      vBot = y;
      iBot = i;}
    if (x > vRight) {
      vRight = x;
      iRight = i;}
    if (x < vLeft) {
      vLeft = x;
      iLeft = i;}
  }
  swap(P[iTop], P[0]); iTop = 0;
  swap(P[iBot], P[1]); iBot = 1;
  swap(P[iLeft], P[2]); iLeft = 2;
  swap(P[iRight], P[3]); iRight = 3;
  auto f0 = newFacetLeft(iTop, P[iLeft], P[iTop]);
  auto f1 = newFacetRight(iTop, P[iTop], P[iRight]);
  auto f2 = newFacetRight(iBot, P[iRight], P[iBot]);
  auto f3 = newFacetLeft(iBot, P[iBot], P[iLeft]);
  f0->next = f1; f1->prev = f0;
  f1->next = f2; f2->prev = f1;
  f2->next = f3; f3->prev = f2;
  f3->next = f0; f0->prev = f3;
  facet* H = f0;

  if(verbose) {
    cout << "initial-hull = ";
    printHull(H, H);
  }

  pointNode* PN = newA(pointNode, n);
  for(intT i=0; i<n; ++i) {
    PN[i] = pointNode(P[i]);
    auto ptr = H;
    do {
      if (ptr->visibleFrom(P[i])) {
        PN[i].seeFacet = ptr;
        ptr->push_back(&PN[i]);
        break; //each point only records one visible facet
      }
      ptr = ptr->next;
    } while (ptr != H);
  }

  //given facet H, find furthest visible point (only among recorded)
  auto findPivot = [&](facet* H)
		   {
		     auto f = H;
		     while (f->size()<=0 && f->next!=H) f = f->next;
		     if (f->size() <= 0) return (pointNode*)NULL;
		     point2d l = f->p1;
		     point2d r = f->p2;
		     auto triangArea = [&](intT idx)
				       {
					 return triArea(l, r, f->at(idx)->p);
				       };
		     intT idx = sequence::maxIndex<double>((intT)0, (intT)f->size(), greater<floatT>(), triangArea);
		     return f->at(idx);
		   };

  intT i = 4;//first four points are in
  while(1) {
    if(verbose) cout << "--- iter" << i << endl;

    //find a pivot point to be processed next
    pointNode pr;
    if (farPivot) {
      pointNode* prp = findPivot(H);
      if (!prp) break;
      pr = *prp;
    } else {
      if (i>=n) break;
      pr = PN[i];
      if (!pr.seeFacet) {
	i++;
	continue;}
    }

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
      facet* new1 = newFacetLeft(i, start->p1, pr.p);
      facet* new2 = newFacetRight(i, pr.p, end->p1);

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
          for(intT j=0; j<ptr->size(); ++j) {
            auto seePt = ptr->at(j);
	    seePt->seeFacet = NULL;
            if (new1->visibleFrom(seePt->p)) {
	      seePt->seeFacet = new1;
              new1->push_back(seePt);
            } else if (new2->visibleFrom(seePt->p)) {
	      seePt->seeFacet = new2;
              new2->push_back(seePt);
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

      if(verbose) {
        cout << "hull = ";
        printHull(H, H);
      }
    }
    i++;
  }

  free(PN);
#ifndef SILENT
  cout << "hull-time = " << t.next() << endl;
#else
  cout << t.next() << endl;
#endif

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

#ifndef SILENT
  intT hSize = 0;
  auto ptr = H;
  do {
    hSize ++;
    ptr = ptr->next;
  } while (ptr != H);
  cout << "hull size = " << hSize << endl;
#endif

  intT* I = newA(intT, n);
  return _seq<intT>(I,1);//dummy
}
