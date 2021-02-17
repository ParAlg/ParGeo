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
#include "pbbs/sampleSort.h"
#include "pbbs/randPerm.h"
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
  pointNode** seeList;
  intT seeListSize;

  // Constructors
  facet(point2d p11, point2d p22): p1(p11), p2(p22) {
    seeList = NULL;
    seeListSize = 0;
  }

  // Methods
  void assign(pointNode** seeIn, intT seeSizeIn) {
    //if (seeSizeIn < 0) abort();
    seeList = seeIn;
    seeListSize = seeSizeIn;}
  bool visibleFrom(point2d p) {return triArea(p1, p2, p) > numericKnob;}
  intT size() {return seeListSize;};
  pointNode* at(intT i) {return seeList[i];}
  void cacheStart(facet* start) {prev = start;}
  void cacheEnd(facet* end) {next = end;}
  facet* getStart() {return prev;}
  facet* getEnd() {return next;}
};

static std::ostream& operator<<(std::ostream& os, const facet f) {
  os << "(" << f.p1 << ", " << f.p2 << ")";
  return os;
}

struct pointNode {
  typedef double floatT;
  point2d p;
  facet* seeFacet;//maintain one edge visible, change from time to time
  pointNode(point2d pp): p(pp), seeFacet(NULL) {};
  pointNode() {};
  inline floatT x() {return p.x();}
  inline floatT y() {return p.y();}
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
  static bool brute = false;//visibility check todo remove
  static bool verify = false;
  static bool farPivot = false;//todo add in farpivot optimization

  auto facets = newA(facet, 2*n);

  auto newFacetLeft = [&](intT i, point2d p11, point2d p22) {
			facets[2*i] = facet(p11, p22);
			return &facets[2*i];
		      };
  auto newFacetRight = [&](intT i, point2d p11, point2d p22) {
			 facets[2*i+1] = facet(p11, p22);
			 return &facets[2*i+1];
		       };
  auto getFacetTmp = [&](intT i) {
			return &facets[2*i];
		      };

  timing t; t.start();
  timing tt; tt.start();

  pointNode* PN = newA(pointNode, n);
  parallel_for(0, n, [&](intT i) {PN[i] = pointNode(P[i]);});

  auto findExtreme = [&] (intT s, intT e, intT& iTop, intT &iBot, intT &iLeft, intT &iRight, floatT& vTop, floatT& vBot, floatT& vLeft, floatT& vRight) {
		       // Akl–Toussaint heuristic
		       vTop = floatMin();
		       vBot = floatMax();
		       vLeft = floatMax();
		       vRight = floatMin();
		       for(intT i=s; i<e; ++i) {
			 floatT x = PN[i].x();
			 floatT y = PN[i].y();
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
		     };

  intT numWorkers = getWorkers();
  intT batchSize = floor(n/numWorkers);
  intT indice[numWorkers*4];
  floatT vals[numWorkers*4];
  parallel_for(0, numWorkers, [&](intT i){
				intT s = i*batchSize;
				intT e;
				if (i < numWorkers-1) e = s+batchSize;
				else e = s+max(batchSize, n-s);
				intT iTopLoc; intT iBotLoc;
				intT iLeftLoc; intT iRightLoc;
				floatT vTopLoc; floatT vBotLoc;
				floatT vLeftLoc; floatT vRightLoc;
				findExtreme(s, e, iTopLoc, iBotLoc, iLeftLoc, iRightLoc,
					    vTopLoc, vBotLoc, vLeftLoc, vRightLoc);
				indice[i*4+0] = iTopLoc; indice[i*4+1] = iBotLoc;
				indice[i*4+2] = iLeftLoc; indice[i*4+3] = iRightLoc;
				vals[i*4+0] = vTopLoc; vals[i*4+1] = vBotLoc;
				vals[i*4+2] = vLeftLoc; vals[i*4+3] = vRightLoc;
			      }, 1);

  intT iTop, iBot, iLeft, iRight;
  floatT vTop = floatMin();
  floatT vBot = floatMax();
  floatT vLeft = floatMax();
  floatT vRight = floatMin();
  for(intT i=0; i<numWorkers; ++i) {
    if (vals[i*4+0]>vTop) iTop = indice[i*4+0];
    if (vals[i*4+1]<vBot) iBot = indice[i*4+1];
    if (vals[i*4+2]<vLeft) iLeft = indice[i*4+2];
    if (vals[i*4+3]>vRight) iRight = indice[i*4+3];
  }

  swap(PN[iTop], PN[0]); iTop = 0;
  swap(PN[iBot], PN[1]); iBot = 1;
  swap(PN[iLeft], PN[2]); iLeft = 2;
  swap(PN[iRight], PN[3]); iRight = 3;
  // Ordering matters
  auto f0 = newFacetLeft(iTop, PN[iLeft].p, PN[iTop].p);
  auto f1 = newFacetRight(iTop, PN[iTop].p, PN[iRight].p);
  auto f2 = newFacetLeft(iBot, PN[iRight].p, PN[iBot].p);
  auto f3 = newFacetRight(iBot, PN[iBot].p, PN[iLeft].p);
  f0->next = f1; f1->prev = f0;
  f1->next = f2; f2->prev = f1;
  f2->next = f3; f3->prev = f2;
  f3->next = f0; f0->prev = f3;
  facet* H = f0;

  if(verbose) {
    cout << "initial-hull = ";
    printHull(H, H);
  }

  // Sort points so those assigned to the same facet are adjacent
  parallel_for(4, n, [&](intT i) {
		       PN[i].seeFacet = NULL;
		       auto ptr = H;
		       do {
			 if (ptr->visibleFrom(PN[i].p)) {
			   PN[i].seeFacet = ptr;
			   break;
			 }
			 ptr = ptr->next;
		       } while (ptr != H);
		     });

  sampleSort(PN+4, n-4, [&](pointNode p1, pointNode p2) {
			  return p1.seeFacet < p2.seeFacet;
			});

  auto SL = newA(pointNode*, n*2);
  parallel_for(0, n, [&](intT i) {SL[i] = &PN[i];});
  parallel_for(n, n*2, [&](intT i) {SL[i] = NULL;});

  // Assign segments of pointers to 4 initial facets
  pair<intT, facet*> assignment[5];
  intT cnt = 0;
  parallel_for(4, n, [&](intT i) {
		       if (i == 4) {
			 if (SL[i]->seeFacet)
			   assignment[utils::fetchAndAdd(&cnt, 1)]
			     = make_pair(i, SL[i]->seeFacet);
		       } else if (SL[i]->seeFacet != SL[i-1]->seeFacet) {
			 assignment[utils::fetchAndAdd(&cnt, 1)]
			   = make_pair(i, SL[i]->seeFacet);
		       }
		     });
  assignment[4] = make_pair(n, (facet*)NULL);
  sort(assignment, assignment+5);

  auto ptr = H;
  do {
    for (int j=0; j<4; ++j) {
      if (assignment[j].second == ptr) {
	ptr->assign(&SL[assignment[j].first],
		    assignment[j+1].first-assignment[j].first);
      }
    }
    ptr = ptr->next;
  } while (ptr != H);

  intT processed = assignment[0].first;//skip points that are already in
  auto pointers = newA(pointNode*, n);
  parallel_for(0, processed, [&](intT i) { pointers[i] = NULL;});
  parallel_for(processed, n, [&](intT i) { pointers[i] = &PN[i];});
  floatT initTime = tt.next();

  randPerm(pointers+processed, n-processed);
  //std::random_shuffle(pointers+processed, pointers+n);
  floatT shuffleTime = tt.stop();

  auto flag = newA(intT, n+1);//packing flag
  auto pointers2 = newA(pointNode*, n);//extra memory for packing

  intT* reservation = newA(intT, 2*n);
  parallel_for(0, 2*n, [&](intT i) {reservation[i]=intMax();});

  pointNode*** seeLists = newA(pointNode**, 2*n);
  parallel_for(0, 2*n, [&](intT i) {seeLists[i]=NULL;});

  auto fIdx = [&](facet* f) {return intT(f-facets);};
  auto pIdx = [&](pointNode* p) {return intT(p-PN);};

  floatT reserveTime = 0;
  floatT confirmTime = 0;
  floatT processTime = 0;
  floatT packTime = 0;

  // Thread private memory for storing new hull
  intT numWorker = num_workers();
  facet** hullStarts = newA(facet*, numWorker);

  floatT batch = 4;

  intT totalRounds = 0;
  intT serialRounds = 0;
  floatT totalTime = 0;
  while(processed < n) {
    timing rt; rt.start();
    totalRounds ++;

    intT b = min(floatT(batch), floatT(n-processed));
    intT s = processed;
    intT roundProcessed;

    if (b < 2000)
      goto serialRound;

    for(intT i=0; i<numWorker; ++i) hullStarts[i] = NULL;

    // Reservation
    tt.start();

    parallel_for(s, s+b, [&](intT i) {
			   pointNode* pr = pointers[i];

			   if (pr->seeFacet) {

			     //find range of visible facets [left, right)
			     pair<facet*, facet*> conflicts;
			     if (brute)
			       conflicts = findVisible(H, pr->p);
			     else
			       conflicts = findVisible(pr->seeFacet, pr->p);
			     facet* start = conflicts.first;
			     facet*   end = conflicts.second;

			     facet* tmp = getFacetTmp(pIdx(pr));
			     tmp->cacheStart(start);
			     tmp->cacheEnd(end);

			     if (start && end) {
			       auto ptr = start->prev;
			       do {
				 utils::writeMin(&reservation[fIdx(ptr)], pIdx(pr));
				 ptr = ptr->next;
			       } while (ptr != end->next);
			     }
			   }
			 });

    reserveTime += tt.next();

    // Confirm reservation

    parallel_for(s, s+b, [&](intT i) {
			   pointNode* pr = pointers[i];

			   if (!pr->seeFacet) {
			     // Note: it's important to memorize the invisibility
			     // of these points here, since during processing
			     // there hull and visibility will be updated
			     pointers[i] = NULL; // Not visible to OLD hull, skip
			   } else {
			     facet* tmp = getFacetTmp(pIdx(pr));
			     facet* start = tmp->getStart();
			     facet* end = tmp->getEnd();

			     if (start && end) {

			       //confirm reservation
			       bool reserved = true;
			       auto ptr = start->prev;
			       do {
				 if (reservation[fIdx(ptr)] == pIdx(pr)) {
				   reservation[fIdx(ptr)] = intMax();//reset
				 } else {
				   reserved = false;
				 }
				 ptr = ptr->next;
			       } while (ptr != end->next);

			       if (reserved) reservation[fIdx(start)] = pIdx(pr);//marker
			     }
			   }
			 });

    confirmTime += tt.next();

    // Processing
    parallel_for(s, s+b, [&](intT i) {

			   pointNode *pr = pointers[i];

			   if (pr) {
			     if(verbose) cout << " pr = " << pr->p << endl;

			     facet* tmp = getFacetTmp(pIdx(pr));
			     facet* start = tmp->getStart();
			     facet* end = tmp->getEnd();
			     if (start && end) {
			       //confirm reservation
			       bool reserved = reservation[fIdx(start)]==pIdx(pr);

			       if (reserved) {
				 reservation[fIdx(start)] = intMax();//late reset
				 //compute new faces
				 //note: the index for new facet memory should be based on PN
				 // since i indexes into pointers array that is reordered
				 // after each round
				 facet* new1 = newFacetLeft(pIdx(pr), start->p1, pr->p);
				 facet* new2 = newFacetRight(pIdx(pr), pr->p, end->p1);

				 pointers[i] = NULL; //will process now, no further actions needed

				 intT cnt = 0;
				 auto ptr = start;
				 do {
				   cnt += ptr->size();
				   ptr = ptr->next;
				 } while (ptr != end);

				 if (cnt > 0) { // Re-assign visible points

				   bool copy = false;
				   ptr = start;
				   pointNode** A = NULL;
				   do {
				     if (!copy && ptr != start) {
				       if (ptr->seeList < ptr->prev->seeList){ //wrapped around
					 copy = true;
					 A = ptr->prev->seeList + ptr->prev->size();
				       }}

				     if (copy) {
				       ptr->seeList = A;
				       if (ptr->size() > 0) {
					 granular_for(0, ptr->size(), 2000,
						      [&](intT j) { A[j] = ptr->at(j);});
					 A += ptr->size();
				       }
				     }
				     ptr = ptr->next;
				   } while (ptr != end);

				   pointNode** SLM = start->seeList;

				   // Update points that see facet* ptr
				   ptr = start;
				   do {
				     granular_for(0, ptr->size(), 2000,
						  [&](intT j) {
						    auto seePt = ptr->at(j);
						    if (new1->visibleFrom(seePt->p)) {
						      seePt->seeFacet = new1;
						    } else if (new2->visibleFrom(seePt->p)) {
						      seePt->seeFacet = new2;
						    } else {
						      seePt->seeFacet = NULL;
						    }
						  });
				     ptr = ptr->next;
				   } while (ptr != end);

				   intT sortCnt;
				   ptr = start;
				   while(ptr != end) {
				     sortCnt = ptr->next->seeList - start->seeList;
				     if (sortCnt < 0)
				       sortCnt = ptr->seeList + ptr->size() - start->seeList;
				     ptr = ptr->next;
				   };

				   // Update new facets' visible points list
				   //  SLM: where the seeLists of affected facets start
				   //  sortCnt: length SLM affected
				   // Now need to split
				   if (true) {
				     sampleSort(SLM, sortCnt, [&](pointNode* p1, pointNode* p2) {
							    return p1->seeFacet < p2->seeFacet;});

				     intT asn[3]; asn[0] = -1; asn[1] = -1;
				     granular_for(0, sortCnt, 2000, [&](intT j) {
								  if (j == 0 || (SLM[j]->seeFacet != SLM[j-1]->seeFacet)) {
								    if (SLM[j]->seeFacet==new1)
								      asn[0] = j;
								    else if (SLM[j]->seeFacet==new2)
								      asn[1] = j;
								  }
								});
				     asn[2] = sortCnt;

				     new1->assign(&SLM[0], 0);
				     if (asn[0]>= 0) {
				       intT new1Size;
				       if (asn[1]>=0) new1Size = asn[1]-asn[0];
				       else new1Size = asn[2]-asn[0];
				       if (new1Size>0)
					 new1->assign(&SLM[asn[0]], new1Size);
				     }

				     if (asn[1]>=0 && asn[2]-asn[1] > 0)
				       new2->assign(&SLM[asn[1]], asn[2]-asn[1]);
				     else
				       new2->assign(new1->seeList + new1->size(), 0);

				   } else {
				     // Two pass split, todo bug
				     bool* F = newA(bool, sortCnt);
				     auto SLM2 = newA(pointNode*, sortCnt);
				     granular_for(0, sortCnt, 2000, [&](intT j) {
								      if (SLM[j]->seeFacet) F[j] = 0;
								      else F[j] = 1; });
				     intT nonEmpty = sequence::split(SLM, SLM2, F, sortCnt);

				     if (nonEmpty > 0) {
				       granular_for(0, nonEmpty, 2000, [&](intT j) {
									 if (SLM2[j]->seeFacet==new1) F[j] = 0;
									 else F[j] = 1;
								       });

				       intT size1 = sequence::split(SLM2, SLM, F, nonEmpty);

				       for(intT j=nonEmpty; j<sortCnt; ++j)
					 SLM[j]->seeFacet = NULL;

				       if (size1>0) new1->assign(&SLM[0], size1);
				       else new1->assign(&SLM[0], 0);
				       if (nonEmpty-size1>0) new2->assign(&SLM[size1], nonEmpty-size1);
				       else new2->assign(new1->seeList+new1->size(), 0);
				     } else {
				       new1->assign(&SLM[0], 0);
				       new2->assign(&SLM[0], 0);
				     }
				     free(SLM2);
				     free(F);
				   }
				 } // End re-assigning visible points

				 if(verbose) cout << "adding = " << *new1 << ", " << *new2 << endl;

				 //update hull
				 start->prev->next = new1; new1->prev = start->prev;
				 new1->next = new2; new2->prev = new1;
				 new2->next = end; end->prev = new2;
				 hullStarts[worker_id()] = new1;

				 if(verbose) {
				   cout << "hull = ";
				   printHull(H, H);
				 }
			       }
			     }
			   }
			 });//end par_for

    for(intT i=0; i<numWorker; ++i) {
      if (hullStarts[i]) {
	H = hullStarts[i];
	break;
      }
    }

    processTime += tt.next();

    //pack unprocessed

    if (b < 2000) {
      intT lPt = s;
      intT rPt = s+b-1;

      while (lPt < rPt) {
        if (pointers[lPt]) {//right
          while (pointers[rPt] && lPt < rPt) {//right
            rPt--;
          }
          if (lPt < rPt) {
            swap(pointers[lPt], pointers[rPt]);
            rPt--; }
          else { break;}
        }
        lPt++;
      }
      if (pointers[lPt]==NULL) lPt++;//left
      roundProcessed = lPt - s;
    } else {
      intT e = s+b;
      parallel_for(s, e, [&](intT i){
			   if (pointers[i]) {
			     if (!pointers[i]->seeFacet) {
			       pointers[i] = NULL;
			       flag[i] = 0;
			     } else {
			       flag[i] = 1;
			     }
			   } else {
			     flag[i] = 0;
			   }
			 });
      flag[e] = sequence::prefixSum(flag, s, e);
      intT numConf = flag[e];
      parallel_for(s, e, [&](intT i) {
			   if (flag[i] != flag[i+1]) {
			     pointers2[flag[i]] = pointers[i];
			   }
			 });
      parallel_for(0, numConf, [&](intT i) {
				 pointers[e-numConf+i] = pointers2[i];
			       });
      roundProcessed = b - numConf;
    }

    packTime += tt.stop();

    // cout << "new processed = " << roundProcessed << "/" << b << ": " << rt.stop() << endl;
    // intT hSize = 0;
    // {
    //   auto ptr = H;
    //   do {
    // 	hSize ++;
    // 	ptr = ptr->next;
    //   } while (ptr != H);
    // }
    totalTime += rt.stop();
    goto finalize;

  serialRound:
    serialRounds ++;
    for(intT i=s; i<s+b; ++i) {

      pointNode *pr = pointers[i];

      if (pr) {
	if (!pr->seeFacet) {
	  pointers[i] = NULL;
	  continue;
	}

	if(verbose) cout << " pr = " << pr->p << endl;

	pair<facet*, facet*> conflicts;
	if (brute)
	  conflicts = findVisible(H, pr->p);
	else
	  conflicts = findVisible(pr->seeFacet, pr->p);
	facet* start = conflicts.first;
	facet*   end = conflicts.second;

	if (start && end) {

	  facet* new1 = newFacetLeft(pIdx(pr), start->p1, pr->p);
	  facet* new2 = newFacetRight(pIdx(pr), pr->p, end->p1);
	  pointers[i] = NULL; //will process now, no further actions needed

	  intT cnt = 0;
	  auto ptr = start;
	  do {
	    cnt += ptr->size();
	    ptr = ptr->next;
	  } while (ptr != end);

	  if (cnt > 0) { // Re-assign visible points
	    bool copy = false;
	    ptr = start;
	    pointNode** A = NULL;
	    do {
	      if (!copy && ptr != start) {
		if (ptr->seeList < ptr->prev->seeList){ //wrapped around
		  copy = true;
		  A = ptr->prev->seeList + ptr->prev->size();
		}}

	      if (copy) {
		ptr->seeList = A;
		if (ptr->size() > 0) {
		  granular_for(0, ptr->size(), 2000,
			       [&](intT j) { A[j] = ptr->at(j);});
		  A += ptr->size();
		}
	      }
	      ptr = ptr->next;
	    } while (ptr != end);

	    pointNode** SLM = start->seeList;

	    // Update points that see facet* ptr
	    ptr = start;
	    do {
	      granular_for(0, ptr->size(), 2000,
			   [&](intT j) {
			     auto seePt = ptr->at(j);
			     if (new1->visibleFrom(seePt->p)) {
			       seePt->seeFacet = new1;
			     } else if (new2->visibleFrom(seePt->p)) {
			       seePt->seeFacet = new2;
			     } else {
			       seePt->seeFacet = NULL;
			     }
			   });
	      ptr = ptr->next;
	    } while (ptr != end);

	    intT sortCnt;
	    ptr = start;
	    while(ptr != end) {
	      sortCnt = ptr->next->seeList - start->seeList;
	      if (sortCnt < 0)
		sortCnt = ptr->seeList + ptr->size() - start->seeList;
	      ptr = ptr->next;
	    };

	    sampleSort(SLM, sortCnt, [&](pointNode* p1, pointNode* p2) {
				       return p1->seeFacet < p2->seeFacet;});

	    intT asn[3]; asn[0] = -1; asn[1] = -1;
	    granular_for(0, sortCnt, 2000, [&](intT j) {
					     if (j == 0 || (SLM[j]->seeFacet != SLM[j-1]->seeFacet)) {
					       if (SLM[j]->seeFacet==new1)
						 asn[0] = j;
					       else if (SLM[j]->seeFacet==new2)
						 asn[1] = j;
					     }
					   });
	    asn[2] = sortCnt;

	    new1->assign(&SLM[0], 0);
	    if (asn[0]>= 0) {
	      intT new1Size;
	      if (asn[1]>=0) new1Size = asn[1]-asn[0];
	      else new1Size = asn[2]-asn[0];
	      if (new1Size>0)
		new1->assign(&SLM[asn[0]], new1Size);
	    }

	    if (asn[1]>=0 && asn[2]-asn[1] > 0)
	      new2->assign(&SLM[asn[1]], asn[2]-asn[1]);
	    else
	      new2->assign(new1->seeList + new1->size(), 0);

	  } // End re-assigning visible points

	  //update hull
	  start->prev->next = new1; new1->prev = start->prev;
	  new1->next = new2; new2->prev = new1;
	  new2->next = end; end->prev = new2;
	  H = new1;
	}
      }
    }//end serial for
    roundProcessed = b;

  finalize:
    processed = s+roundProcessed;
    intT numConflicts = batch - roundProcessed;

    if (numConflicts < 0.5 * batch) {
      batch *= 1.2;
    } else {
      batch /= 10;
    }
    if (batch <= 0) batch = 10;

  }//end while

  parallel_for(0, 2*n, [&](intT i) {
			 // Freeing see lists allocated
			 if (seeLists[i]) free(seeLists[i]);
		       });
  free(SL);
  free(seeLists);
  free(PN);
  free(hullStarts);
  free(pointers2);
  free(flag);

#ifndef SILENT
  cout << "rounds = " << totalRounds << endl;
  cout << " serial-rounds = " << serialRounds << endl;
  cout << "hull-time = " << t.next() << endl;
  cout << " shuffle-time = " << shuffleTime << endl;
  cout << " init-time = " << initTime << endl;
  cout << " reserve-time = " << reserveTime << endl;
  cout << " confirm-time = " << confirmTime << endl;
  cout << " process-time = " << processTime << endl;
  cout << " pack-time = " << packTime << endl;
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
  {
    intT hSize = 0;
    auto ptr = H;
    do {
      hSize ++;
      ptr = ptr->next;
    } while (ptr != H);
    cout << "hull size = " << hSize << endl;
  }
#endif

  free(facets);
  free(reservation);
  free(pointers);

  intT* I = newA(intT, n);
  return _seq<intT>(I,1);//dummy
}
