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

namespace randInternal {

  struct pointNode;

  struct facet {
    // Data fields
    facet* next;
    facet* prev;
    intT p1;
    intT p2;//(p1->p2) is clockwise
    vector<pointNode*>* seeList;//points that can see facet

    // Constructors
    facet(intT p11, intT p22): p1(p11), p2(p22) {
      seeList = new vector<pointNode*>();}//todo allocation make more efficient (might call in parallel)

    facet(): p1(-1), p2(-1), seeList(NULL), next(NULL), prev(NULL) {}

    // Methods
    bool visibleFrom(point2d* P, intT p) {return triArea(P[p1], P[p2], P[p]) > numericKnob;}
    void push_back(pointNode* p) {seeList->push_back(p);}
    intT size() {return seeList->size();};
    pointNode* at(intT i) {return seeList->at(i);}
  };

  static std::ostream& operator<<(std::ostream& os, const facet f) {
    os << "(" << f.p1 << ", " << f.p2 << ")";
    return os;
  }

  struct pointNode {
    intT p;
    facet* seeFacet;//maintain one edge visible, change from time to time
    //don't have to maintain all facets, can search pretty easily
    pointNode(intT pp): p(pp), seeFacet(NULL) {};
    pointNode() {};
  };

  pair<facet*, facet*> findVisible(point2d* P, facet* head, intT p) {
    while (head->prev->visibleFrom(P, p)) head = head->prev;

    auto ptr = head;
    facet* start = NULL;
    intT n = 0;
    do {
      if (n <= 0) {
	if (ptr->visibleFrom(P, p)) {
	  n++;
	  start = ptr;
	}
      } else {
	if (!ptr->visibleFrom(P, p)) {
	  return make_pair(start, ptr);
	}
	n++;
      }
      ptr = ptr->next;
    } while (ptr != head);
    return make_pair((facet*)NULL, (facet*)NULL);
  }
}

using namespace randInternal;

// Internal rand hull, takes in a partial hull H
// - array facets, if allocated needs to be at least of size 2*n
facet* randHullInternalSerial(point2d* P, pointNode* PN, intT n, facet* H, facet* facets=NULL) {
  static bool verbose = false;
  static bool farPivot = true;
  static bool pivotSample = true;
  static intT sampleSize = 1000;

  if (!facets) {
    auto facets = newA(facet, 2*n);}

  auto newFacetLeft = [&](intT i, intT p11, intT p22) {
			facets[2*i] = facet(p11, p22);
			return &facets[2*i];
		      };
  auto newFacetRight = [&](intT i, intT p11, intT p22) {
			 facets[2*i+1] = facet(p11, p22);
			 return &facets[2*i+1];
		       };

  timing t; t.start();

  //given facet H, find furthest visible point (only among recorded)
  auto findPivot = [&](facet* H)
		   {
		     auto f = H;
		     while (f->size()<=0 && f->next!=H) f = f->next;
		     if (f->size() <= 0) return (pointNode*)NULL;
		     point2d l = P[f->p1];
		     point2d r = P[f->p2];
		     auto triangArea = [&](intT idx)
				       {
					 return triArea(l, r, P[f->at(idx)->p]);
				       };
		     intT idx = sequence::maxIndex<double>((intT)0, (intT)f->size(), greater<floatT>(), triangArea);
		     return f->at(idx);
		   };
  auto findPivotSample = [&](facet* H)
			 {
			   auto f = H;
			   while (f->size()<=0 && f->next!=H) f = f->next;
			   if (f->size() <= 0) return (pointNode*)NULL;
			   point2d l = P[f->p1];
			   point2d r = P[f->p2];
			   auto triangArea = [&](intT idx)
					     {
					       return triArea(l, r, P[f->at(idx)->p]);
					     };
			   intT reductionSize = min((intT)f->size(), sampleSize);
			   intT idx = sequence::maxIndex<double>((intT)0, reductionSize, greater<floatT>(), triangArea);
			   return f->at(idx);
			 };

  intT i = 0;
  intT roundCount = 0;
  while(1) {
    if(verbose) cout << "--- iter" << i << endl;

    //find a pivot point to be processed next
    pointNode pr;
    if (farPivot) {
      if (pivotSample) {
	pointNode* prp = findPivotSample(H);
	if (!prp) break;
	pr = *prp;
      } else {
	pointNode* prp = findPivot(H);
	if (!prp) break;
	pr = *prp;
      }
    } else {
      if (i>=n) break;
      pr = PN[i];
      if (!pr.seeFacet) {
	i++;
	continue;}
    }

    roundCount ++;

    if(verbose) cout << " pr = " << pr.p << endl;

    //find range of visible facets [left, right)
    pair<facet*, facet*> conflicts;
    conflicts = findVisible(P, pr.seeFacet, pr.p);

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

      // update visibility pointers
      auto ptr = start;
      do {
	//update points that see them
	for(intT j=0; j<ptr->size(); ++j) {
	  auto seePt = ptr->at(j);
	  seePt->seeFacet = NULL;
	  if (new1->visibleFrom(P, seePt->p)) {
	    seePt->seeFacet = new1;
	    new1->push_back(seePt);
	  } else if (new2->visibleFrom(P, seePt->p)) {
	    seePt->seeFacet = new2;
	    new2->push_back(seePt);//todo can i just push the index
	  }
	}
	ptr = ptr->next;
      } while (ptr != end);

      if(verbose) cout << "adding = " << *new1 << ", " << *new2 << endl;

      //update hull
      start->prev->next = new1; new1->prev = start->prev;
      new1->next = new2; new2->prev = new1;
      new2->next = end; end->prev = new2;
      H = new1;

    }
    i++;
  }
  if(verbose) cout << "#rounds = " << roundCount << endl;

  return H;
}

// Start with an externally computed partial-hull
// - P[I[0,n)] stores points to be processed
// - P[Itmp[0,m)] stores a partial convex hull
// - Resulting hull will be stored in I (overwrite)
// - Hull size will be returned
intT randHullExternalSerial(point2d* P, intT* I, intT n, intT* Itmp, intT m) {
  facet* H = newA(facet, m+2*n);
  facet* facets = H+m;

  for(intT i=0; i<m; ++i) {
    H[i].p1 = Itmp[i];
    H[i].p2 = Itmp[(i+1)%m];
    H[i].seeList = new vector<pointNode*>();//todo allocation
    H[i].prev = &H[(i+m-1)%m];
    H[i].next = &H[(i+1)%m];
  }

  pointNode* PN = newA(pointNode, n);
  for(intT i=0; i<n; ++i) {
    PN[i] = pointNode(I[i]);
    auto ptr = H;
    do {
      if (ptr->visibleFrom(P, I[i])) {
        PN[i].seeFacet = ptr;
        ptr->push_back(&PN[i]);
        break; //each point only records one visible facet
      }
      ptr = ptr->next;
    } while (ptr != H);
  }

  auto newH = randHullInternalSerial(P, PN, n, H, facets);

  intT mm = 0;
  auto head = newH;
  do {
    if (head->p1 == Itmp[0]) break;
    head = head->next;
  } while (head != newH);

  auto ptr = head;
  do {
    if (ptr - H < m || ptr->p1 == Itmp[0]) {
      ptr = ptr->next;
    } else {
      I[mm++] = ptr->p1;
      ptr = ptr->next;
    }
  } while (ptr != head);

  free(H);
  free(PN);
  return mm;
}

// Call from a partial hull
// I contains indices for the existing hull, it is required that
// it has enough space to hold the entire hull (size n)
// I will be overwritten
intT randHullSerial(point2d* P, intT n, intT* I, intT m) {
  auto H = newA(facet, m+2*n);
  for(intT i=0; i<m; ++i) {
    H[i].p1 = I[i];
    H[i].p2 = I[(i+1)%m];
    H[i].seeList = new vector<pointNode*>();//todo allocation
    H[i].prev = &H[(i+m-1)%m];
    H[i].next = &H[(i+1)%m];
  }

  auto facets = H+m;

  pointNode* PN = newA(pointNode, n);//todo consider allocate outside
  for(intT i=0; i<n; ++i) {
    PN[i] = pointNode(i);
    auto ptr = H;
    do {
      if (ptr->visibleFrom(P, i)) {
        PN[i].seeFacet = ptr;
        ptr->push_back(&PN[i]);
        break; //each point only records one visible facet
      }
      ptr = ptr->next;
    } while (ptr != H);
  }

  auto newH = randHullInternalSerial(P, PN, n, H, facets);

  intT mm = 0;
  auto ptr = newH;
  do {
    I[mm++] = ptr->p1;
    ptr = ptr->next;
  } while (ptr != newH);
  free(H);

  return mm;
}

// Start from scratch, initialize with 4 extremum
intT randHullSerial(point2d* P, intT n, intT* I) {
  // Akl–Toussaint heuristic
  // picks the left, right, top, bottom points to form an initial hull of size 4

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

  I[0] = iLeft;
  I[1] = iTop;
  I[2] = iRight;
  I[3] = iBot;

  return randHullSerial(P, n, I, 4);
}

// Start from scratch, with a random facet
intT randHullNaiveSerial(point2d* P, intT n, intT* I) {
  I[0] = 0;
  I[1] = 1;

  return randHullSerial(P, n, I, 2);
}
