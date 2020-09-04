// This code is part of the project "A Parallel Batch-Dynamic Data Structure
// for the Closest Pair Problem"
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

#ifndef SIEVE_GRID_H
#define SIEVE_GRID_H

#include <math.h>
#include <limits>
#include <vector>
#include "geometry.h"
#include "kdTree.h"
#include "shared.h"
#include "pbbs/sequence.h"
#include "pbbs/ndHash.h"
#include "pbbs/sampleSort.h"
#include "pbbs/quickSort.h"
using namespace std;

template<int dim> struct cellHash;

/**
 *  A cell class that represents each box in the grid.
 */
template<int dim>
struct cell {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef cell<dim> cellT;

  pointT *P;//point array of the cell (stored externally)
  pointT coordP;//average coordinate of the cell
  intT numPoints;

  /**
  *   Constructor 1, that initializes an empty cell.
  */
  cell(): P(NULL), numPoints(-1) {}

  /**
  *   Constructor 2, that only initializes the coordinate of the cell base on an external point.
  *   @param PP external point.
  */
  cell(pointT PP): coordP(PP), numPoints(1) {}

  /**
  *   Constructor 3, that initializes the cell with an external point array.
  *   @param PP external point array.
  *   @param nn size of PP.
  */
  cell(pointT *PP, intT nn): P(PP), numPoints(nn) {}

  /**
  *   Given pMin, grid size r, computes the cell coordinate (center).
  *   @param pMin global minimum coordinate (all dimensions).
  *   @param r grid cell side length.
  */
  void computeCoord(pointT pMin, double r) {
    for(int i=0; i<dim; ++i) {
      coordP.updateCoordinate(i, r/2+pMin[i]+(double)floor((P[0][i]-pMin[i])/r)*r);
    }
  }

  void init() {P=NULL; numPoints=-1;}

  inline bool isEmpty() {return numPoints<0;}

  intT size() {return numPoints;}

  floatT *coordinate() {
    if (isEmpty()) return NULL;
    else return coordP.coordinate();}

  floatT coordinate(int i) {return coordP.coordinate(i);}

  pointT operator[](intT i) {
    if (i >= size()) {
      cout << "error, cell access out of bound, abort" << endl; abort();}
    return P[i];}

  /**
  *   Determines if another cell is bordering on the current.
  *   @param c2 another cell.
  *   @param r cell side length.
  *   @return boolean.
  */
  bool nextTo(cellT& c2, floatT r) {
    for (int i=0; i<dim; ++i) {
      if (abs(coordinate(i) - c2.coordinate(i)) > r*1.000001) {
        return false;
      }
    }
    return true;}
};

/**
 *  A grid class that stores points in axis-aligned box cells.
 */
template<int dim>
struct grid {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef pointPair<dim> pointPairT;
  typedef cell<dim> cellT;
  typedef Table<cellHash<dim>,intT> tableT;

  static const bool noRandom = false;
  static const unsigned int prime = -5;
  static const unsigned int mask = -1;
  static const unsigned int range = (1 << 29);

  floatT r;
  pointT pMin;
  cellT *cells;
  intT numCells, n;
  tableT* table;
  kdTree<dim, cellT>* tree = NULL;

  int rands[10] = {846930886, 1681692777, 1714636915, 1957747793, 424238335, 719885386, 1649760492, 596516649, 1189641421, 120120309};
  int randInt[dim];
  uintT modprime(intT* x, intT n) {
    unsigned long long temp = 0;
    uintT key = 0;
    for (intT i=0; i<n; i++) {
      temp = (long long) x[i] * (long long) randInt[i];
      temp = (temp & mask) + 5 * (temp >> 32);
      if (temp >= prime) temp -= prime;
      temp += key;
      if (temp >= prime) temp -= prime;
      key = (uintT) temp;
    }
    return key;
  }

  uintT hash(floatT *x) {
    intT xx[dim];
    for(int i=0; i<dim; ++i) {
      xx[i] = (intT) floor((x[i]-pMin[i])/r);}
    return modprime(xx, dim);}

  inline int compareCoordinate(floatT* x1, floatT* x2) {
    for(int i=0; i<dim; ++i) {
      intT xx1 = (intT) floor((x1[i]-pMin[i])/r);
      intT xx2 = (intT) floor((x2[i]-pMin[i])/r);
      if (xx1 != xx2) {
        if (xx1 > xx2) return 1;
        else return -1;
      }}
    return 0;
  }

  inline void resetCells() {
    par_for(intT i=0; i<n; ++i) {
      cells[i].init();
    }
    numCells = 0;
  }

  inline void initRand() {
    srand(time(NULL));
    for (intT i = 0; i < dim; i++) {
      if(noRandom) randInt[i] = rands[i] % range + 1;
      else randInt[i] = rand() % range + 1;}
  }

  /**
   *  Constructor.
   *  @param nn maximum allowed number of points.
   *  @param pMinn global coordinate minimum.
   *  @param rr grid size length.
   *  Constructor.
   */
  grid(intT nn, pointT pMinn, floatT rr): pMin(pMinn), r(rr), n(nn) {
    initRand();
    cells = newA(cellT, n);
    resetCells();
  }

  ~grid() {free(cells);}

  void addTable(tableT* tablee) {table = tablee;}

  void updateR(floatT rr) {r = rr;}

  /**
   *  Iterates the cell neighborhood of the cell given coordinate.
   *  @param center center coordinate (coordP) of the cell.
   *  @param doTerm a function that takes in cell pointer, processes it, and return whether to stop the procedure.
   */
  template<class func>
  void cellNeighborhood(cellT* center, func& doTerm) {
    if (center->isEmpty()) {
      cout << "error, finding neighbor for empty cell, abort"<< endl; abort();}
    intT numNeighbor = 3;//3**dim
    for(int i=1; i<dim; ++i) numNeighbor*=3;
    pointT p = pointT();
    cellT bait = cellT(p);
    for(int i=0; i<numNeighbor; ++i) {
      int ii = i;
      for(int j=0; j<dim; ++j) {
        bait.coordP.x[j] = center->P[0][j]+((floatT)(ii%3)-1)*r*1.0000001;//numerical stability
        ii /= 3;}
      cellT* nbr = table->find(&bait);
      if (!nbr->isEmpty()) {
        if (doTerm(nbr)) break;}
    }
  }

  /**
   *  Serially, clears the grid, and re-insert new points.
   *  @param P new point array.
   *  @param nn size of P.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return number of cells created.
   */
  intT reInsertSerial(pointT* P, intT nn, bool useTree=false) {
    resetCells();
    auto pLess = [&] (pointT a, pointT b) {
                   return pointGridCmp<dim, pointT>(a, b, pMin, r);};
    quickSortSerial(P, nn, pLess);

    numCells = 0;
    cells[numCells].P = &P[0];
    for (intT i=1; i<nn; ++i) {
      intT cellSize = 0;
      if (compareCoordinate(P[i].x, P[i-1].x)) {
        cells[numCells].numPoints = &P[i] - cells[numCells].P;
        cells[numCells].computeCoord(pMin, r);
        numCells++;
        cells[numCells].P = &P[i];}
    }
    cells[numCells].numPoints = &P[nn] - cells[numCells].P;
    cells[numCells].computeCoord(pMin, r);
    numCells++;

    if (!useTree) {
      table->setActive(min(n,numCells*3));
      table->clear();
      for (intT i=0; i<numCells; ++i) {
        table->insert(&cells[i]);}
    }
    tree = NULL;
    return numCells;
  }

  /**
   *  In parallel, clears the grid, and re-insert new points.
   *  @param P new point array.
   *  @param nn size of P.
   *  @param flag auxiliary memory of size nn.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return number of cells created.
   */
  intT reInsertParallel(pointT* P, intT nn, intT* flag=NULL, bool useTree=false) {
    bool freeFlag=false;
    if(!flag) {
      flag=newA(intT,nn+1);
      freeFlag=true;}
    resetCells();
    auto pLess = [&] (pointT a, pointT b) {
                   return pointGridCmp<dim, pointT>(a, b, pMin, r);};
    sampleSort(P, nn, pLess);
    flag[0] = 1;
    par_for(intT i=1; i<nn; ++i) {
      if (compareCoordinate(P[i].x, P[i-1].x)) flag[i] = 1;
      else flag[i] = 0;
    }
    numCells = sequence::prefixSum(flag, 0, nn);

    flag[nn] = numCells;
    cells[0].P = P;
    par_for(intT i=1; i<nn; ++i) {
      intT ii = i+1;
      if (flag[ii] != flag[ii-1]) {
        cells[flag[ii]-1].P = &P[ii-1];}
    }
    par_for(intT i=0; i<numCells-1; ++i) {
      cells[i].numPoints = cells[i+1].P - cells[i].P;
      cells[i].computeCoord(pMin, r);
    }
    cells[numCells-1].numPoints = &P[nn] - cells[numCells-1].P;
    cells[numCells-1].computeCoord(pMin, r);

    if (!useTree) {
      table->setActive(min(n,numCells*3));
      table->clear();
      par_for(intT i=0; i<numCells; ++i) {
        table->insert(&cells[i]);
      }
    }
    tree = NULL;
    if(freeFlag) free(flag);
    return numCells;
  }

  inline bool isSparse(pointT pp) {
    cellT* target = table->find(new cellT(pp));
    if (target->isEmpty()) {
      cout << "error, isSparse empty target, abort" << endl; abort();}

    intT numPoints = 0;
    auto doTerm = [&](cellT *c) {
                    numPoints += c->size();
                    return !(numPoints <= 1);};

    cellNeighborhood(target, doTerm);
    return numPoints <= 1;
  }

  inline bool isSparseKdTree(pointT pp) {
    pointT p = pointT(pp);
    cellT bait = cellT(&p,1);
    bait.computeCoord(pMin, r);

    intT numPoints = 0;
    auto term = [&]() {return numPoints > 1;};//prunes tree traversal
    auto doTerm = [&](cellT *c) {//prunes leaf traversal
                    if (c->nextTo(bait, r)) {
                      numPoints += c->size();
                      return numPoints > 1;
                    }};

    tree->rangeNeighbor(&bait, r*sqrt(dim), term, doTerm);//check range
    return numPoints <= 1;
  }

  /**
   *  Serially, collect the dense points from the newly inserted points.
   *  @param P points newly inserted just now.
   *  @param PP auxiliary array of size nn.
   *  @param nn size of P.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return sequence built on PP, containing the dense points.
   */
  _seq<pointT> checkInsertSerial(pointT* P, pointT* PP, intT nn, bool useTree=false) {
    if (useTree) tree = new kdTree<dim, cellT>(cells, numCells, false);
    intT ii = 0;
    for (intT i=0; i<nn; ++i) {
      if (useTree) {
        if (!isSparseKdTree(P[i])) PP[ii++] = P[i];
      } else {
        if (!isSparse(P[i])) PP[ii++] = P[i];}
    }
    return _seq<pointT>(PP,ii);
  }

  /**
   *  In parallel, collect the dense points from the newly inserted points.
   *  @param P points newly inserted just now.
   *  @param PP auxiliary array of size nn.
   *  @param nn size of P.
   *  @param flag auxiliary memory of size nn.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return sequence built on PP, containing the dense points.
   */
  _seq<pointT> checkInsertParallel(pointT* P, pointT* PP, intT nn, bool* flag, bool useTree=false) {
    if (useTree) tree = new kdTree<dim, cellT>(cells, numCells, true);
    par_for(intT i=0; i<nn; ++i) {
      if (useTree) {
        if (isSparseKdTree(P[i])) flag[i] = false;
        else flag[i] = true;
      } else {
        if (isSparse(P[i])) flag[i] = false;
        else flag[i] = true;}
    }
    auto pCreate = [&](intT i) {return pointT(P[i].x);};
    auto seqPP = sequence::pack(PP, flag, (intT)0, nn, pCreate);
    return seqPP;
  }

  /**
   *  Serially, computes the closest pair in the grid.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return closest pair.
   */
  pointPairT closestPairSerial(bool useTree=false) {
    if (useTree && tree == NULL) tree = new kdTree<dim, cellT>(cells, numCells, false);
    auto R = pointPairT();
    for (intT i=0; i<numCells; ++i) {
      auto c = &cells[i];
      vector<pointT> P;
      auto gather = [&](cellT *c) {
                      for (intT i=0; i<c->size(); ++i) {
                        P.push_back(c->P[i]);}
                      return false;};
      if (useTree) {
        tree->rangeNeighbor(c, r*sqrt(2), [&](){return false;}, gather);
      } else {
        cellNeighborhood(c, gather);
      }
      auto myCp = pointPairT();
      for (intT i=0; i<P.size(); ++i) {
        for (intT j=i+1; j<P.size(); ++j) myCp.closer(P[i],P[j]);
      }
      R.closer(myCp);
    }
    return R;
  }

  /**
   *  In parallel, computes the closest pair in the grid.
   *  @param useTree whether uses a kdTree to help neighborhood queries.
   *  @return closest pair.
   */
  pointPairT closestPairParallel(bool useTree=false) {
    if (useTree && tree == NULL) tree = new kdTree<dim, cellT>(cells, numCells, true);
    auto Rs = newA(pointPairT, numCells);
    par_for(intT i=0; i<numCells; ++i) {
      auto c = &cells[i];
      vector<pointT> P;
      auto gather = [&](cellT *c) {
                      for (intT i=0; i<c->size(); ++i) {
                        P.push_back(c->P[i]);}
                      return false;};
      if (useTree) {
        auto notComplete = [&]() {return false;};
        tree->rangeNeighbor(c, r*sqrt(2), [&](){return false;}, gather);
      } else {
        cellNeighborhood(c, gather);}
      auto myCp = pointPairT();
      for (intT i=0; i<P.size(); ++i) {
        for (intT j=i+1; j<P.size(); ++j) myCp.closer(P[i],P[j]);
      }
      Rs[i] = myCp;
    }
    auto getDist = [&](intT i){return Rs[i].dist;};
    intT iMin = sequence::minIndex<floatT,intT>(0, numCells, getDist);
    auto R = Rs[iMin]; free(Rs);
    return R;
  }
};

template<int dim>
struct cellHash {
  typedef point<dim> pointT;
  typedef grid<dim> gridT;
  typedef cell<dim>* eType;
  typedef cell<dim>* kType;
  gridT* grids;
  cellHash(gridT* gridss): grids(gridss) {};
  eType empty() {return new cell<dim>();}
  kType getKey(eType v) {return v;}
  uintT hash(kType c) {
    // if (c->isEmpty()) {
    //   cout << "error, hashing empty cell, abort" << endl; abort();}
    return grids->hash(c->coordinate());
  }
  int cmp(kType c1, kType c2) {
    if (c1->isEmpty() || c2->isEmpty()) return 1;
    return grids->compareCoordinate(c1->coordinate(), c2->coordinate());
  }
  bool replaceQ(eType c1, eType c2) {return 0;}
};

#endif
