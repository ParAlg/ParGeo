// This code is part of the project "Fast Parallel Algorithms for
// Euclidean Minimum SpanningTree and Hierarchical Spatial Clustering"
// Copyright (c) 2020 Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
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

#ifndef DENDROGRAM_H
#define DENDROGRAM_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "pbbs/parallel.h"
#include "pbbs/utils.h"
#include "pbbs/sampleSort.h"
#include "pbbs/gettime.h"

using namespace std;

//#define CHECKER
//#define VERBOSE
//#define DEBUG_VERBOSE

namespace dendrogram {
  struct edge {
    intT p, q;
    double weight;};

  struct edgePointer {
    intT p, idx;};

  struct eulerNode {
    eulerNode* oppo;
    inline eulerNode* getPrev(eulerNode* tour, intT m) {
      return &tour[(m+idx-1)%m];
    }
    inline eulerNode* getNext(eulerNode* tour, intT m) {
      return &tour[(idx+1)%m];
    }

    intT from, to, rank, marked;
    intT idx;
    intT tmpFrom, tmpTo;
    intT dis = 0;
    inline void setNext(intT next) {
      tmpFrom = next;}
    inline intT next() {
      return tmpFrom;}

    void printNode() {
      cout << "(" << from << "," << to << ")";}

    intT l=-1;
    inline void updateLabel(intT ll) {l=ll;};
    inline bool smallerLabel(intT ll) {
      if (ll<l) {
        l=ll;
        return true;
      } else {
        false;
      }
    }
    inline intT label() {return l;};
    inline intT hasLabel() {return l>=0;};
    intT head=-1;
    inline void setHead(intT spSize) {head=spSize;};
    inline void unsetHead() {head=-1;};
    inline intT isHead() {return head>=1;};
    inline intT spSize() {return head;};

    bool isEqual(eulerNode* eN) {
      if (eN->from == from && eN->to == to) return true;
      return false;}

    bool isFromAbut(eulerNode* eN) {
      if (eN->from == from || eN->to == from) {
        return !isEqual(eN);
      }
      return false;
    }
    bool isToAbut(eulerNode* eN) {
      if (eN->from == to || eN->to == to) {
        return !isEqual(eN);
      }
      return false;
    }

    void printTour(eulerNode* tour, intT m) {
      tour[0].printNode(); cout << ": ";
      for (intT i=1; i<m; ++i) {
        tour[i].printNode();cout << " ";}
      cout << "| ";
      tour[0].printNode();
    }

    vector<eulerNode*> fromNbrs(eulerNode* tour, intT m) {
      vector<eulerNode*> nbrs;
      eulerNode* eN = getPrev(tour,m);
      while (!eN->isEqual(oppo)) {
        if (isFromAbut(eN)) {
          nbrs.push_back(eN);
          eN = eN->getPrev(tour,m);
        } else {
          eN = eN->oppo->getPrev(tour,m);}
      }
      return nbrs;
    }

    vector<eulerNode*> toNbrs(eulerNode* tour, intT m) {
      vector<eulerNode*> nbrs;
      eulerNode* eN = getNext(tour,m);
      while (!eN->isEqual(oppo)) {
        if (isToAbut(eN)) {
          nbrs.push_back(eN);
          eN = eN->getNext(tour,m);
        } else {
          eN = eN->oppo->getNext(tour,m);
        }
      }
      return nbrs;
    }

    eulerNode* isLead(intT* vtxDis, eulerNode* tour, intT m) {
      if (vtxDis[from] > vtxDis[to]) return NULL;
      auto fromAbut = fromNbrs(tour, m);
      eulerNode* eN;
      bool hasPrev = false;
      for (intT i=0; i<fromAbut.size(); ++i) {
        eN = fromAbut[i];
        if (eN->to == from && vtxDis[eN->from] < vtxDis[eN->to]) {
          hasPrev = true;
          if (eN->marked) return eN;}
      }
      if (hasPrev) return NULL;//has unqualified light prev edge
      else return oppo;//light edge has start node
    }

    void propLabel(eulerNode* stopNode, eulerNode* tour, intT m) {
      intT linkVertex = to;
      eulerNode* eN = getNext(tour, m);
      intT count = 0;
      while (true) {
        if (!eN->marked && linkVertex == eN->from) {
          eN->updateLabel(label());
          count++;
          linkVertex = eN->to;}
        if (eN == stopNode) break;
        else eN = eN->getNext(tour, m);
      }
    }

    void initTmpTour() {
      tmpFrom=-1;
      tmpTo=-1;}

    void tmpTour() {
      if(tmpFrom>=0) from=tmpFrom;
      if(tmpTo>=0) to=tmpTo;}

    void changeLabel(intT newId, eulerNode* tour, intT m) {
      eulerNode* eN;
      auto fromAbut = fromNbrs(tour,m);
      for (intT i=0; i<fromAbut.size(); ++i) {
        eN = fromAbut[i];
        if (eN->marked && eN->to == from) eN->tmpTo = newId;
      }
      auto toAbut = toNbrs(tour,m);
      for (intT i=0; i<toAbut.size(); ++i) {
        eN = toAbut[i];
        if (eN->marked && eN->from == to) eN->tmpFrom = newId;}
    }
    double weight;
  };

  struct eulerWeightCmp {
    bool operator() (eulerNode i, eulerNode j) {
      return i.weight < j.weight;}
  };

  struct spPackIdxCmp {
    intT* idx;
    eulerNode* tour;
  spPackIdxCmp(intT* idxx, eulerNode* tourr, intT m):idx(idxx), tour(tourr) {
    par_for(intT i=0; i<m; ++i) idx[i]=i;}
    bool operator() (intT ii, intT jj) {//input: index in tour
      eulerNode i = tour[ii];
      eulerNode j = tour[jj];
      if (i.label() < 0 && j.label() >= 0) return true;
      if (i.label() >= 0 && j.label() < 0) return false;
      if (i.label() < 0 && j.label() < 0) return i.rank < j.rank;
      if (i.label() == j.label()) return i.rank<j.rank;
      return i.label() < j.label();}
  };

  struct eulerRankIdxCmp {
    intT* idx;
    eulerNode* tour;
  eulerRankIdxCmp(intT* idxx, eulerNode* tourr, intT m):idx(idxx), tour(tourr) {
    par_for(intT i=0; i<m; ++i) idx[i]=i;}
    bool operator() (intT ii, intT jj) {//input: index in tour
      eulerNode i = tour[ii];
      eulerNode j = tour[jj];
      return i.rank < j.rank;}
  };

  struct edgePointerCmp {
    bool operator() (const edgePointer i, const edgePointer j) {
      return (i.p < j.p);}
  };

  struct ufNode {
    intT father, size, ptr;};

  struct unionFind {
    intT n;
    vector<ufNode> A;
  unionFind(intT nn): n(nn) {
    A.resize(n);
    for (intT i=0; i<n; i++) {
      A[i].father = i;
      A[i].ptr = -1;
      A[i].size = 1;}
  }
    inline intT father(intT p) {
      return A[p].father == p ? p : A[p].father = father(A[p].father);}
  };

  struct dendroNode {
  dendroNode(): weight(-1){};
    void setEmpty() {weight=-1;}
    bool isEmpty() {return weight<0;}
    double weight;//edge weight
    intT vtx[2];//edge vertices
    intT ptr[2];};//offset (in ddg) to children

  // the main class
  template<class edgeT>
    struct directedDendro{
      static const intT startNode = 0;
      intT n;//incoming vertex count
      vector<edge> edges;//input
      vector<edgePointer> edgePointers;//auxiliary
      vector<intT> vtxPtr;//auxiliary
      vector<intT> vtxDis;//stores distance to source mst vertex
      vector<eulerNode> eulerTour;
      vector<dendroNode> ddg;//output

      //constructor and main function
    directedDendro(edgeT *edgesIn, intT nn): n(nn) {
      //printMstCsv(edges, n-1);
      //init
      edgePointers.resize(2*n-2);
      eulerTour.resize(2*n-2);
      vtxPtr.resize(n+1);
      vtxDis.resize(n);
      ddg.resize(n-1);
      edges.resize(2*n-2);
      generateTree(edgesIn, n);

      // construct eulerTour
      semisort();
      constructEulerTour();
      listRanking();

      // construct ddg
#ifdef CHECKER
      par_for(intT i=0; i<n-1; ++i) ddg[i].setEmpty();
#endif

      // construct ddg: init
      unionFind uf = unionFind(n);
      eulerNode* tour = &eulerTour[0];
      intT m = eulerTour.size();

      eulerNode* newTour = newA(eulerNode, m);
      intT* newOrder = newA(intT, m*2+1);
      intT* newOrderInverse = newOrder+m;

      // construct ddg: preprocessing, sort euler tour such that linked edges are adjacent
      // basically sort by eulerNode.rank, but need to maintain oppo pointer, so is more complicated
      // oppo pointer: basically (a,b).oppo = (b,a)
      par_for(intT i=0; i<m; ++i) tour[i].idx = i;
      sampleSort(newOrder, m, eulerRankIdxCmp(newOrder, tour, m));
      par_for(intT i=0; i<m; ++i) {
        newTour[i] = tour[newOrder[i]];
        newOrderInverse[newOrder[i]] = i;}//where tour[i] went in newTour
      par_for(intT i=0; i<m; ++i) {
        newTour[i].oppo = &newTour[newOrderInverse[tour[newOrder[i]].oppo->idx]];}
      swap(tour, newTour);

      double* weights = newA(double, m);

      // construct ddg: main function
      constructDDGParallel(tour, m, uf, &ddg[0], 0, weights, newOrder, newTour);
      //constructDDGSerial(tour, m, uf, &ddg[0], 0);

#ifdef CHECKER
      for(intT i=0; i<n-1; ++i) {
        if (ddg[i].isEmpty()) {
          cout << "loc = " << i << endl;
          cout << "error, checking ddg, found empty, abort" << endl;abort();
        }
      }
      // printDDG();
#endif
    }

      inline void addEdge(intT i, intT x, intT y, double w) {
        edge e;
        e.p = x; e.q = y; e.weight = w;
	edges[i] = e;
        e.p = y; e.q = x;
	edges[i+1] = e;
      }

      inline void generateTree(edgeT *edgess, intT nn) {
        n = nn;
        for (intT i = 0; i < n-1; i++) {
          auto myEdge = edgess[i];
          addEdge(i*2, myEdge.u, myEdge.v, myEdge.weight);
	}
      }

      inline void semisort() {
        par_for (intT i = 0; i < edges.size(); i++) {
          edgePointers[i].p = edges[i].p;
          edgePointers[i].idx = i;
          eulerTour[i].from = edges[i].p;
          eulerTour[i].to = edges[i].q;
          eulerTour[i].weight = edges[i].weight;
        }
        sampleSort(&edgePointers.at(0), edgePointers.size(), edgePointerCmp());
        par_for (intT i = 0; i < edges.size(); i++) {
          eulerTour[edgePointers[i].idx].idx = i;}
        edgePointer temp;
        temp.p = n; temp.idx = -1;
        edgePointers.push_back(temp);
      }

      inline void constructEulerTour() {
        vtxPtr[0] = 0;
        par_for (intT i = 0; i < 2 * n - 2; i++) {
          if (edgePointers[i].p != edgePointers[i + 1].p) {
            vtxPtr[edgePointers[i + 1].p] = i + 1;}
        }
        par_for (intT i = 0; i < n; i++) {
          for (intT j = vtxPtr[i]; j < vtxPtr[i + 1]; j++) {
            intT last = j;
            if (last == vtxPtr[i]) last = vtxPtr[i + 1];
            last--;
            intT ptr = edgePointers[j].idx;
            eulerTour[edgePointers[last].idx ^ 1].setNext(ptr);
            eulerTour[ptr].oppo = &eulerTour[ptr ^ 1];
          }
        }
      }

      void listRanking() {
        int curEdge = 0, curDis = 0;
        const int SEGMAX = 2000;
        int sum[SEGMAX], mark[SEGMAX], next[SEGMAX];
        int SEG = SEGMAX;
        if (SEG > n / 1000) SEG = n / 1000;
        if (SEG == 0) SEG = 1;
        int startPos = edgePointers[vtxPtr[startNode]].idx;
        mark[0] = startPos;
        eulerTour[startPos].marked = 1;
        for (int lp = 1; lp < SEG; lp++) {
          do
            mark[lp] = (long long)lp * 462387463 % (2 * n - 2);
          while (eulerTour[mark[lp]].marked);
          eulerTour[mark[lp]].marked = lp + 1;
        }
        par_for (int lp = 0; lp < SEG; lp++) {
          if (eulerTour[mark[lp]].marked == lp + 1) {
            int curRank = 1, curPtr = mark[lp];
            eulerTour[curPtr].rank = 0;
            curPtr = eulerTour[curPtr].next();
            while (eulerTour[curPtr].marked == 0) {
              eulerTour[curPtr].rank = curRank;
              curRank++;
              curPtr = eulerTour[curPtr].next();
            }
            sum[lp] = curRank;
            next[lp] = eulerTour[curPtr].marked - 1;
          }
          else sum[lp] = 0;
        }
        int curSum = sum[0];
        sum[0] = 0;
        int curPtr = next[0];
        while (curPtr != 0) {
          int temp = sum[curPtr];
          sum[curPtr] = curSum;
          curSum += temp;
          curPtr = next[curPtr];
        }
        par_for (int lp = 0; lp < SEG; lp++) {
          if (eulerTour[mark[lp]].marked == lp + 1) {
            int curPtr = mark[lp];
            eulerTour[curPtr].rank = sum[lp];
            curPtr = eulerTour[curPtr].next();
            while (eulerTour[curPtr].marked == 0) {
              eulerTour[curPtr].rank += sum[lp];
              curPtr = eulerTour[curPtr].next();
            }
          }
        }

        par_for (int i = 0; i < n - 1; i++) {
          if (eulerTour[i * 2].rank < eulerTour[i * 2 + 1].rank) {
            eulerTour[i * 2].dis = 1;
            eulerTour[i * 2 + 1].dis = -1;
          }
          else {
            eulerTour[i * 2 + 1].dis = 1;
            eulerTour[i * 2].dis = -1;
          }
        }

        par_for (int lp = 0; lp < SEG; lp++) {
          if (eulerTour[mark[lp]].marked == lp + 1) {
            int curDis = eulerTour[mark[lp]].dis, curPtr = mark[lp];
            curPtr = eulerTour[curPtr].next();
            while (eulerTour[curPtr].marked == 0) {
              eulerTour[curPtr].dis += curDis;
              curDis = eulerTour[curPtr].dis;
              curPtr = eulerTour[curPtr].next();
            }
            sum[lp] = curDis;
          }
          else sum[lp] = 0;
        }
        curSum = sum[0];
        sum[0] = 0;
        curPtr = next[0];
        while (curPtr != 0) {
          int temp = sum[curPtr];
          sum[curPtr] = curSum;
          curSum += temp;
          curPtr = next[curPtr];
        }
        par_for (int lp = 0; lp < SEG; lp++) {
          if (eulerTour[mark[lp]].marked == lp + 1) {
            int curPtr = mark[lp];
            eulerTour[curPtr].dis += sum[lp];
            vtxDis[eulerTour[curPtr].to] = eulerTour[curPtr].dis;
            curPtr = eulerTour[curPtr].next();
            while (eulerTour[curPtr].marked == 0) {
              eulerTour[curPtr].dis += sum[lp];
              vtxDis[eulerTour[curPtr].to] = eulerTour[curPtr].dis;
              curPtr = eulerTour[curPtr].next();
            }
          }
        }

      }

      void constructDDGSerial(eulerNode* tour, intT m, unionFind& uf, dendroNode* dendro, intT dendroAbs) {
#ifdef VERBOSE
        cout << "\n  constructDDGSerial -------- dendro[" << dendroAbs << "]-?" << endl;
        cout << "  tour size = " << m << endl;
        cout << "  [input tour] ";tour->printTour(tour, m);cout << endl;
#endif

        if (m%2>0) {
          cout << "  error, odd length tour = "<< m<<", abort" << endl;abort();}

#ifdef CHECKER
        for (intT i=0; i<m/2; ++i) {
          if (!dendro[i].isEmpty()) {
            cout << "  loc = " << i + dendroAbs << endl;
            cout << "  error, ddg basecase overwrites ddg, abort" << endl;abort();
          } else {
            dendro[i].weight = 1;}//set non-empty
        }
#endif

        quickSortSerial(tour, m, eulerWeightCmp());
        intT ii = 0;
        for (intT i=0; i<m; i++) {
          if (tour[i].from < tour[i].to) continue;//process one side
          intT p = tour[i].from;
          intT q = tour[i].to;
          if (vtxDis[tour[i].from] > vtxDis[tour[i].to]) swap(p, q);
          dendro[ii].vtx[0] = p;
          dendro[ii].vtx[1] = q;
          dendro[ii].weight = tour[i].weight;
          dendro[ii].ptr[0] = uf.A[uf.father(p)].ptr;
          dendro[ii].ptr[1] = uf.A[uf.father(q)].ptr;
          uf.A[uf.father(p)].ptr = dendroAbs+i;
          uf.A[uf.father(q)].father = uf.father(p);
          uf.A[uf.father(p)].size += uf.A[uf.father(q)].size;
          ii++;
        }
      }

      void constructDDGParallel(eulerNode* tour, intT m, unionFind& uf, dendroNode* dendro,
                                intT dendroAbs, double* weights, intT* newOrder, eulerNode* newTour) {
#ifdef VERBOSE
        cout << "\nconstructDDGParallel -------- dendro[" << dendroAbs << "]-?" << endl;
        cout << "tour size = " << m << endl;
        cout << "[input tour] ";tour->printTour(tour, m);cout << endl;
#endif

        if (m%2>0) {
          cout << "error, odd length tour = "<< m<<", abort" << endl;abort();}

        if (m<=(n*0.5)) {
          constructDDGSerial(tour, m, uf, dendro, dendroAbs);
          return;
        }
#ifdef DEBUG_VERBOSE
        printPartialMstCsv(&tour[0], m);
#endif

        par_for(intT i=0; i<m; ++i) {
          weights[i] = tour[i].weight;}

        sampleSort(weights, m, less<double>());
        intT numLight = 0.9*m;
        while (numLight && weights[numLight-1]==weights[numLight]) numLight++;
        if (numLight%2>0) numLight++;

        double wMedian = weights[numLight];
        intT numHeavy = m-numLight;
#ifdef VERBOSE
        cout << "weightMedian = " << wMedian << endl;
        cout << "numHeavy = " << numHeavy << endl;
        cout << "numLight = " << numLight << endl;
#endif

        //mark heavy edge as 1
        par_for(intT i=0; i<m; ++i) {
          eulerNode* eN = &tour[i];
          eN->unsetHead();
          eN->updateLabel(-1);
          eN->idx = i;
          if (eN->weight >= wMedian) {
            eN->marked = 1;
          } else {
            eN->marked = 0;
          }
        }

        //identify light subproblems
        //if identified lead edge of subproblem, propagate the label to all its light edges
        par_for(intT i=0; i<m; ++i) {
          eulerNode* eN = &tour[i];
          if (!eN->marked) {
            eulerNode* leadHeavy = eN->isLead(&vtxDis[0], tour, m);
            if (leadHeavy != NULL) {
              eN->updateLabel(eN->from);
              eN->propLabel(leadHeavy->oppo, tour, m);}
          }
        }

        //reconnect heavy edges into one subproblem
        //replace removed edge-endpoints to lead vertex of removed subproblem
        par_for(intT i=0; i<m; ++i) {
          eulerNode* eN = &tour[i];
          if (eN->marked) eN->initTmpTour();
        }
        par_for(intT i=0; i<m; ++i) {
          eulerNode* eN = &tour[i];
          if (!eN->marked) {
            eN->changeLabel(eN->label(), tour, m);
          }
        }
        par_for(intT i=0; i<m; ++i) {
          eulerNode* eN = &tour[i];
          if (eN->marked) eN->tmpTour();
        }

        // packing light subproblems
        intT* newOrderInverse = newOrder+m;
        sampleSort(newOrder, m, spPackIdxCmp(newOrder, tour, m));
        par_for(intT i=0; i<m; ++i) {
          newTour[i] = tour[newOrder[i]];
          newOrderInverse[newOrder[i]] = i;}//where tour[i] went in newTour
        par_for(intT i=0; i<m; ++i) {
          newTour[i].oppo = &newTour[newOrderInverse[tour[newOrder[i]].oppo->idx]];}
        swap(tour, newTour);

        //determine subproblem sizes (spSizes)
#ifdef DEBUG_VERBOSE
        for(intT i=0; i<m; ++i) {
          if (tour[i].label()>=0) continue;
          tour[i].printNode();cout << " id: " << tour[i].label() << " ";
          if (tour[i].marked) cout << " mark ";
          cout << endl;
        }
#endif
        par_for(intT i=1; i<m; ++i) {
          if (tour[i].label()>=0) {
#ifdef DEBUG_VERBOSE
            tour[i].printNode();cout << " id: " << tour[i].label() << " ";
#endif
            if (tour[i-1].label()!=tour[i].label()) {
              intT spSize = 1;
              while(i+spSize<m && tour[i+spSize].label()==tour[i].label()) {
                tour[i+spSize].unsetHead();
                spSize++;}
              tour[i].setHead(spSize);
#ifdef DEBUG_VERBOSE
              cout  << "x";
              cout << " size: " << spSize << " ";
#endif
            }
#ifdef DEBUG_VERBOSE
            cout << endl;
#endif
          }
        }

#ifdef VERBOSE
        cout << "[heavy problem] size = " << numHeavy << " ";
        tour[0].printTour(tour, numHeavy);
        cout << endl;
        for(intT i=0; i<m; ++i) {
          if (tour[i].isHead()) {
            cout << "[light problem] size = " << tour[i].spSize() << " ";
            tour[i].printTour(tour+i, tour[i].spSize());
            cout << endl;}
        }
#endif

        //determine offsets (in *tour) of subproblems
        intT* flag = newOrder;
        flag[0]=1;
        par_for(intT i=1; i<m; ++i) {
          if (tour[i].isHead()) {
            flag[i]=1;
          } else {
            flag[i]=0;
          }
        }
        intT numProb = sequence::prefixSum(flag, 0, m);
        intT* probOffset= newA(intT, numProb);
        par_for(intT i=0; i<m-1; ++i) {
          if (flag[i] != flag[i+1]) probOffset[flag[i]] = i;}
        if (flag[m-1] != numProb) {
          probOffset[flag[m-1]] = m-1;}

        //recursive calls to constructDDG
        par_for(intT i=0; i<numProb; ++i) {
          if (i==0) {//heavy problem
            constructDDGParallel(tour, numHeavy, uf, dendro,
                                 dendroAbs, weights, newOrder, newTour);
          } else {//light problem(s)
            eulerNode* head = tour+probOffset[i];
            constructDDGParallel(head, head->spSize(), uf, dendro+probOffset[i]/2,
                                 dendroAbs+probOffset[i]/2, weights+probOffset[i], newOrder+probOffset[i]*2,
                                 newTour+probOffset[i]);
          }
        }
        free(probOffset);
      }

      void printVtxPtr() {
        cout << "[vtxPtr] ";
        for (intT i=0; i<vtxPtr.size(); ++i) {
          cout << vtxPtr[i] << "  ";}
        cout << endl;
      }

      void printVtxDis() {
        cout << "[vtxDis] ";
        for (intT i=0; i<vtxDis.size(); ++i) {
          cout << i<<":"<<vtxDis[i] << "  ";}
        cout << endl;
      }

      void printEdges() {
        cout << "[edges] ";
        for (intT i=0; i<edges.size(); ++i) {
          cout << "<" << edges.at(i).p << "," << edges.at(i).q
               << ">" << ":" << edges.at(i).weight << "  ";}
        cout << endl;
      }

      void printEdgePointers() {
        cout << "[edgePointers] ";
        for (intT i=0; i<edgePointers.size(); ++i) {
          cout << "<" << edgePointers.at(i).p << ">" << ":" << edgePointers.at(i).idx << "  ";}
        cout << endl;
      }

      void printDDG() {
        cout << "[ddg] <v1,v2>:weight";
        for (int i = 0; i < ddg.size(); ++ i) {
          cout << "<" << ddg[i].vtx[0] << "," << ddg[i].vtx[1] << ">:"
               << ddg[i].weight << "  ";}
        cout << endl;
      }

      void printMstCsv(edgeT* mst, intT m) {
        cout << "mst csv ----- start" << endl;
        for (intT i=0; i<m; ++i) {
          cout << mst[i].u << "," << mst[i].v << "," << mst[i].weight << endl;}
        cout << endl;
        cout << "mst csv ----- end" << endl;
      }

      void printPartialMstCsv(eulerNode* start, intT m) {
        cout << "mst csv ----- start" << endl;
        if (start->from < start->to) {
          cout << start->from << "," << start->to << "," << start->weight << endl;
        }
        eulerNode* eN = start->getNext(start,m);
        while (!start->isEqual(eN)) {
          if (eN->from < eN->to) {
            cout << eN->from << "," << eN->to << "," << eN->weight;
            cout<<endl;
          }
          eN = eN->getNext(start,m);
        }
        cout << "mst csv ----- end" << endl;
      }

    };

}//namespace


#endif
