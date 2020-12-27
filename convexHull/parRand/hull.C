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
#include "ringBuf.h"
using namespace std;

struct pointNode;

struct facet {
  facet* next;
  facet* prev;
  point2d p1;
  point2d p2;//(p1->p2) is clockwise
  ringBuffer<pointNode*>* seeList;
  intT s, e;
  inline intT size() {return e-s;}
  void push_back(pointNode* pn) {
    seeList->at(e++) = pn;}
  pointNode* operator[](intT i) {
    return seeList->at(s+i);
  }
  pointNode*& at(intT i) {
    return seeList->at(s+i);
  }
  facet() {}
  facet(point2d p11, point2d p22, ringBuffer<pointNode*> *sl): p1(p11), p2(p22), seeList(sl) {}
  void init(point2d p11, point2d p22, ringBuffer<pointNode*> *sl) {
    p1 = p11;
    p2 = p22;
    seeList = sl;
  }
  void copy(facet* f) {
    p1 = f->p1;
    p2 = f->p2;
    seeList = f->seeList;
    s = f->s;
    e = f->e;
  }
  bool visibleFrom(point2d p) {return triArea(p1, p2, p) > 1e-9;}//todo numerical stability
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

void printHull(facet* start, facet* end, facet* H=NULL, intT n=-1) {
  cout << "PRINT HULL" << endl;
  auto ptr = start;
  do {
    cout << ptr << ":" << *ptr << " ";
    if (ptr->next->prev != ptr) {
      cout << "next->prev link error" << endl;
      cout << " next = " << ptr->next << endl;
      cout << " next->prev = " << ptr->next->prev << endl;
      abort();
    }
    if (ptr->prev->next != ptr) {
      cout << "prev->next link error" << endl;
      cout << " prev = " << ptr->prev << endl;
      cout << " prev->next = " << ptr->prev->next << endl;
      abort();
    }
    if (H) {
      if (ptr-H >= n) {
	cout << "error, access out of hull bound: " << ptr-H << endl;
	abort();
      }
    }
    ptr = ptr->next;
  } while (ptr != end);
  cout << endl;
}

facet* hullIncInPlace(facet* start, facet* end, facet* new1, facet*new2, facet* hullMem, intT listN) {
  intT lmSize = 0;
  pointNode** lm = &start->at(0);

  auto ptr = start;
  do {
    for(intT i=0; i<ptr->size(); ++i) {
      auto pn = ptr->at(i);
      pn->seeFacet = NULL;
      if (new1->visibleFrom(pn->p) ||
	  new2->visibleFrom(pn->p)) {
	lm[lmSize++] = pn;}
    }
    ptr = ptr->next;
  } while (ptr != end);

  auto predicate = [&](intT i) {
		     return new1->visibleFrom(lm[i]->p);// && !new2->visibleFrom(lm[i]->p);
		   };

  intT lPt = 0;
  intT rPt = lmSize-1;
  while (lPt < rPt) {
    if (predicate(lPt)) {
      lm[lPt]->seeFacet = new1;
    } else {
      while (!predicate(rPt) && lPt < rPt) {
	lm[rPt]->seeFacet = new2;
	rPt--;
      }
      lm[lPt]->seeFacet = new2;
      if (lPt < rPt) {
	lm[rPt]->seeFacet = new1;
	swap(lm[lPt], lm[rPt]);
	rPt--; }
      else {break;}
    }
    lPt++;
  }
  if (predicate(lPt)) {
    lm[lPt]->seeFacet = new1;
    lPt++;
  } else {
    lm[lPt]->seeFacet = new2;
  }

  new1->s = start->s;
  new1->e = new1->s + lPt;
  new2->s = new1->e;
  new2->e = new2->s + lmSize-lPt;

  auto prev = start->prev;
  start->copy(new1);
  start->prev = prev;
  prev->next = start;

  facet* new2new = start + start->size() + 1;
  intT offset = new2new - hullMem;
  if (offset >= listN)//wrap around
    new2new = hullMem + (offset % listN);
  new2new->copy(new2);

  start->next = new2new;
  new2new->prev = start;
  new2new->next = end;
  end->prev = new2new;

  //todo simplify this part
  for(intT i=0; i<start->size(); ++i)
    start->at(i)->seeFacet = start;
  for(intT i=0; i<new2new->size(); ++i)
    new2new->at(i)->seeFacet = new2new;

  return start;
}

facet* hullIncrement(facet* start, facet* end, facet* new1, facet*new2, facet* hullMem, intT listN) {
  intT lmSize = 0;
  auto ptr = start;
  do {
    lmSize += ptr->size();
    ptr = ptr->next;
  } while (ptr != end);

  auto lm2 = newA(pointNode*, lmSize);
  for(intT i=0; i<lmSize; ++i) {
    lm2[i] = NULL;
  }
  intT lPt = 0;
  intT rPt = lmSize-1;

  ptr = start;
  do {
    for(intT j=0; j<ptr->size(); ++j) {
      auto seePt = ptr->at(j);
      seePt->seeFacet = NULL;

      if (new1->visibleFrom(seePt->p)) {
	seePt->seeFacet = new1;
	lm2[lPt++] = seePt;
      } else if (new2->visibleFrom(seePt->p)) {
	seePt->seeFacet = new2;
	lm2[rPt--] = seePt;
      }
    }
    ptr = ptr->next;
  } while (ptr != end);

  new1->s = start->s;
  new1->e = new1->s + lPt;
  for(intT i=0; i<lPt; ++i) {
    new1->at(i) = lm2[i];
  }
  new2->s = start->s + rPt+1;
  new2->e = new2->s + lmSize-rPt-1;
  for(intT i=0; i<lmSize-rPt-1; ++i) {
    new2->at(i) = lm2[i+rPt+1];
  }
  free(lm2);

  auto prev = start->prev;
  start->copy(new1);
  start->prev = prev;
  prev->next = start;

  facet* new2new = start + start->size() + 1;
  intT offset = new2new - hullMem;
  if (offset >= listN)//wrap around
    new2new = hullMem + (offset % listN);
  new2new->copy(new2);

  start->next = new2new;
  new2new->prev = start;
  new2new->next = end;
  end->prev = new2new;

  //todo simplify this part
  for(intT i=0; i<start->size(); ++i)
    start->at(i)->seeFacet = start;
  for(intT i=0; i<new2new->size(); ++i)
    new2new->at(i)->seeFacet = new2new;

  return start;
}

//given facet H, find furthest visible point (only among recorded)
inline pointNode* findPivot(facet* H) {
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
}

facet* incrementParallel(pointNode* PN, intT n, facet* H, facet* hullMem, ringBuffer<pointNode*>& listMem, intT listN) {
  static bool brute = false;//bruteforce takes O(h) time to determine visible facets
  static bool verbose = false;
  static bool farPivot = false;
  auto idx = [&](facet* f) {return f-hullMem;};

  intT* reservation = newA(intT, listN);
  parallel_for(0, listN, [&](intT i) {reservation[i]=intMax();});

  intT* flag = newA(intT, n+1);

  pointNode** pointers = newA(pointNode*, n);
  pointNode** pointers2 = newA(pointNode*, n);
  parallel_for(0, n, [&](intT i){pointers[i] = &PN[i];});

  while(n > 10000) {

    //reservation
    parallel_for(0, n,
		 [&](intT i) {
		   auto pr = *pointers[i];
		   if (pr.seeFacet) {//skip non-conflicting points

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
		       auto ptr = start->prev;
		       do {
			 utils::writeMin(&reservation[idx(ptr)], i);
			 ptr = ptr->next;
		       } while (ptr != end->next);
		     }
		   };
		 });

    //processing
    parallel_for(0, n,
		 [&](intT i) {
		   auto pr = *pointers[i];
		   if (pr.seeFacet) {//skip non-conflicting points
		     if(verbose) cout << " pr = " << pr.p << ", sees " << pr.seeFacet << endl;

		     //find range of visible facets [left, right), AGAIN. todo make efficient
		     pair<facet*, facet*> conflicts;
		     if (brute) {
		       conflicts = findVisible(H, pr.p);
		     } else {
		       conflicts = findVisible(pr.seeFacet, pr.p);
		     }
		     facet* start = conflicts.first;
		     facet*   end = conflicts.second;

		     //confirm reservation
		     bool reserved = true;
		     auto ptr = start->prev;
		     do {
		       if (reservation[idx(ptr)] != i) {
			 reserved = false;
			 break;}

		       ptr = ptr->next;
		     } while (ptr != end->next);

		     if (start && end) {
		       if(reserved) {
			 //compute new faces
			 facet new1 = facet(start->p1, pr.p, &listMem);
			 facet new2 = facet(pr.p, end->p1, &listMem);

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
			   auto tmp = hullIncrement(start, end, &new1, &new2, hullMem, listN);
			   if (i == 0) H = tmp;
			   reservation[idx(H)] = intMax();
			   reservation[idx(H->next)] = intMax();
			   //todo inplace hull inc
			 }

			 flag[i] = 0;//processed
		       } else {
			 flag[i] = 1;//reservation fail
		       }
		     } else {
		       flag[i] = 0;//inside hull
		     }
		   } else {
		     flag[i] = 0;
		   }
		 });//end second for

    flag[n] = sequence::prefixSum(flag, 0, n);
    cout << "processed " << n-flag[n] << "/" << n << endl;
    parallel_for(0, n, [&](intT i){
			 if(flag[i] != flag[i+1]) {
			   pointers2[flag[i]] = pointers[i];
			 }
		       });
    swap(pointers, pointers2);

    n = flag[n];
  }//end while

  free(reservation);
  free(pointers2);
  free(flag);

  //processes remaining serially
  intT i=0;
  cout << "process-serial = " << n << endl;
  while(i < n) {
    auto pr = *pointers[i];
    i++;

    if (pr.seeFacet) {//skip non-conflicting points

      //find range of visible facets [left, right), AGAIN. todo make efficient
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
	facet new1 = facet(start->p1, pr.p, &listMem);
	facet new2 = facet(pr.p, end->p1, &listMem);

	if (!brute) {
	  H = hullIncrement(start, end, &new1, &new2, hullMem, listN);
	  //todo inplace hull inc
	}
      }
    }
  }//end while

  free(pointers);
  return H;
}

facet* incrementSerial(pointNode* PN, intT n, facet* H, facet* hullMem, ringBuffer<pointNode*>& listMem, intT listN) {
  static bool brute = false;//bruteforce takes O(h) time to determine visible facets
  static bool verbose = false;
  static bool localPivot = false;

  intT i = 0;
  while(1) {
    if(verbose) cout << "--- iter" << i << endl;

    //find a pivot point to be processed next
    pointNode pr;
    if (localPivot) {
      pointNode* prp = findPivot(H);
      if (!prp) break;
      pr = *prp;
    } else {//random pivot
      //cout << "i check = " << i << ", n = " << n << endl;
      if (i>=n) break;
      pr = PN[i];
      if (!pr.seeFacet) {
	i++;
	continue;}
    }
    i++;

    if(verbose) cout << " pr = " << pr.p << ", sees " << pr.seeFacet << endl;

    //find range of visible facets [left, right)
    pair<facet*, facet*> conflicts;
    if (brute) {
      conflicts = findVisible(H, pr.p);
    } else {
      conflicts = findVisible(pr.seeFacet, pr.p);
    }

    facet* start = conflicts.first;
    facet*   end = conflicts.second;
    //need to reserve [start, end)
    if (start && end) {
      //compute new faces
      facet new1 = facet(start->p1, pr.p, &listMem);
      facet new2 = facet(pr.p, end->p1, &listMem);

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
        //H = hullIncrement(start, end, &new1, &new2, hullMem, listN);
        H = hullIncInPlace(start, end, &new1, &new2, hullMem, listN);
      }

      if(verbose) {
        cout << "hull = ";
        printHull(H, H);
      }
    }
  }

  return H;
}

_seq<intT> hull(point2d* P, intT n) {
  static bool verbose = false;
  static bool verify = true;

  timing t; t.start();

  auto p0 = P[0].x() < P[1].x() ? P[0] : P[1];//left of p1
  auto p1 = P[0].x() < P[1].x() ? P[1] : P[0];
  auto p2 = P[2];

  ringBuffer<pointNode*> listMem = ringBuffer<pointNode*>(n*3);
  for(intT i=0; i<n*3; ++i) listMem[i] = NULL;

  //todo start with a half-constructed hull
  facet f0, f1, f2;
  facet* H;
  if (triArea(p0, p1, p2) > 0.0) {
    f0.init(p0, p2, &listMem);
    f1.init(p2, p1, &listMem);
    f2.init(p1, p0, &listMem);
    f0.next = &f1; f1.prev = &f0;
    f1.next = &f2; f2.prev = &f1;
    f2.next = &f0; f0.prev = &f2;
    H = &f0;
  } else {
    f0.init(p0, p1, &listMem);
    f1.init(p1, p2, &listMem);
    f2.init(p2, p0, &listMem);
    f0.next = &f1; f1.prev = &f0;
    f1.next = &f2; f2.prev = &f1;
    f2.next = &f0; f0.prev = &f2;
    H = &f0;
  }

  if(verbose) {
    cout << "initial-hull = ";
    printHull(H, H);
  }

  H->s = 0;
  H->e = 0;
  H->next->s = n;
  H->next->e = n;
  H->next->next->s = 2*n;
  H->next->next->e = 2*n;

  pointNode* PN = newA(pointNode, n-3);
  for(intT i=0; i<n-3; ++i) {
    PN[i] = pointNode(P[i+3]);
    auto ptr = H;
    do {
      if (ptr->visibleFrom(P[i+3])) {
        PN[i].seeFacet = ptr;
        ptr->push_back(&PN[i]);
        break; //each point only records one visible facet
      }
      ptr = ptr->next;
    } while (ptr != H);
  }

  facet* H2 = newA(facet, n);
  intT o1 = 1+H->size();
  intT o2 = o1+1+H->next->size();
  H2[0].copy(H);
  H2[0].prev = &H2[o2];
  H2[0].next = &H2[o1];
  H2[o1].copy(H->next);
  H2[o1].prev = &H2[0];
  H2[o1].next = &H2[o2];
  H2[o2].copy(H->next->next);
  H2[o2].prev = &H2[o1];
  H2[o2].next = &H2[0];

  granular_for(0, H2[0].size(), 2000,
	       [&](intT i) {H2[0][i]->seeFacet = &H2[0];});
  granular_for(0, H2[o1].size(), 2000,
	       [&](intT i){H2[o1][i]->seeFacet = &H2[o1];});
  granular_for(0, H2[o2].size(), 2000,
	       [&](intT i){H2[o2][i]->seeFacet = &H2[o2];});

  H = H2;

  if(verbose) {
    auto ptr = H;
    do {
      cout << "facet size = " << ptr->size() << endl;
      ptr = ptr->next;
    } while (ptr != H);
    cout << "init complete" << endl;
  }
  cout << "init-time = " << t.next() << endl;

  //H = incrementSerial(PN, n-3, H, H, listMem, n);
  H = incrementParallel(PN, n-3, H, H, listMem, n);
  free(PN);

  cout << "increment-time = " << t.next() << endl;

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
