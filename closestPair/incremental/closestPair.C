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

#include <math.h>
#include <limits>
#include <vector>
#include "geometry.h"
#include "kdTree.h"
#include "shared.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/ndHash.h"
#include "pbbs/sampleSort.h"
#include "pbbs/quickSort.h"
#include "pbbs/gettime.h"
#include "pbbs/randPerm.h"
using namespace std;

// *************************************************************
//   Indexed point
// *************************************************************

template<int dim>
struct iPoint {
  typedef double floatT;
  typedef point<dim> pointT;
  intT i;
  pointT p;
  iPoint(pointT pp, intT ii): p(pp), i(ii) {}
  iPoint(pointT pp): p(pp), i(-1) {}
  iPoint(): i(-1) {}
  bool isEmpty() {return i<0;}
  floatT operator[](intT i) {return p[i];}
  floatT pointDist(iPoint q) {return p.pointDist(q.p);}
  floatT pointDist(pointT q) {return p.pointDist(q);}
  intT idx() {return i;}
  intT idx(intT ii) {i=ii;}
  pointT pt() {return p;}
  floatT* x() {return p.x;}
  floatT x(int i) {return p.x[i];}
  void x(int i, floatT val) {p.x[i]=val;}
};

template <int dim>
static std::ostream& operator<<(std::ostream& os, const iPoint<dim> v) {
  os << "(" << v.i << ") " << v.p; return os;}

// *************************************************************
//   Parallel hashtable to store grids
// *************************************************************

template<int dim> struct cellHash;

/**
 *  A cell class that represents each box in the grid.
 */
template<int dim>
struct cell {
  typedef double floatT;
  typedef iPoint<dim> iPointT;
  typedef point<dim> pointT;
  typedef cell<dim> cellT;
  typedef Table<cellHash<dim>,intT> tableT;

  iPointT *P;
  pointT coordP;
  intT numPoints;

  cell(): P(NULL), numPoints(-1) {}
  cell(iPointT *PP, intT nn): P(PP), numPoints(nn) {}
  cell(pointT PP): coordP(PP), numPoints(1) {}

  void computeCoord(pointT pMin, double r) {
    for(int i=0; i<dim; ++i) {
      coordP.x[i] = r/2+pMin[i]+(double)floor((P[0][i]-pMin[i])/r)*r;}
  }

  void init() {P=NULL; numPoints=-1;}

  inline bool isEmpty() {return numPoints<0;}

  intT size() {return numPoints;}

  floatT *coordinate() {
    if (isEmpty()) return NULL;
    else return coordP.x;}

  floatT coordinate(int i) {return coordP.x[i];}

  intT idx() {return P[0].idx();}

  iPointT operator[](intT i) {
    if (i >= size()) {
      cout << "error, cell access out of bound, abort" << endl; abort();}
    return P[i];}
};

template<int dim>
struct conflictFinder {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef cell<dim> cellT;
  pointT query;
  floatT r;
  intT conflicts;
  intT threshold;
  conflictFinder(pointT queryy, floatT rr)
  :query(queryy), r(rr), conflicts(0), threshold(3) {
    for (int i=1; i<dim; ++i) threshold *= 3;}
  bool hasConflict() {return conflicts > 1;};
  bool checkComplete(cellT *c) {
    if (c->size() > threshold) {
      conflicts = 2;
      return true;}
    for (intT i=0; i<c->size(); ++i) {
      floatT dist = query.pointDist((*c)[i].pt());
      if (dist == 0) {
        conflicts++;
        if (conflicts > 1) return true;
      } else if (dist < r) {
        conflicts = 2;
        return true;}
    }
    return false;}
  bool isComplete() {return conflicts > 1;}
};

template<int dim>
struct closestFinder {
  typedef double floatT;
  typedef point<dim> iPointT;
  typedef pointPair<dim> pointPairT;
  typedef cell<dim> cellT;
  iPointT query;
  intT duplicates;
  intT threshold;
  pointPairT R;
  closestFinder(iPointT queryy):query(queryy),duplicates(0),threshold(3) {
    for (int i=1; i<dim; ++i) threshold *= 3;
    R = pointPairT();}
  pointPairT getClosest() {return R;}
  bool checkComplete(cellT *c) {
    if (c->size() > threshold) {
      for (intT i=0; i<c->size(); ++i) {
        for (intT j=i+1; j<c->size(); ++j)
          R.closer((*c)[i].pt(), (*c)[j].pt());
      }
      return true;}
    for (intT i=0; i<c->size(); ++i) {
      floatT dist = query.pointDist((*c)[i].pt());
      if (dist == 0) {
        duplicates++;
        if (duplicates > 1) {
          R.closer(query,(*c)[i].pt(),dist);
          return true;}
      } else R.closer(query,(*c)[i].pt(),dist);
    }
    return false;}
  bool isComplete() {return duplicates>1;}
};

/**
 *  A grid class that stores points in axis-aligned box cells.
 */
template<int dim>
struct grid {
  typedef double floatT;
  typedef point<dim> pointT;
  typedef iPoint<dim> iPointT;
  typedef pointPair<dim> pointPairT;
  typedef cell<dim> cellT;
  typedef Table<cellHash<dim>,intT> tableT;
  static const int gridSizeRatio = 1;
  static const bool noRandom = false;
  pointPairT R;
  pointT pMin;
  cellT *cells;
  intT numCells, n;
  tableT* table;

  static const unsigned int prime = -5;
  static const unsigned int mask = -1;
  static const unsigned int range = (1 << 29);
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
    return key;}

  uintT hash(floatT *x) {
    intT xx[dim];
    for(int i=0; i<dim; ++i) {
      xx[i] = (intT) floor((x[i]-pMin[i])/r());}
    return modprime(xx, dim);}

  inline void resetCells() {
    par_for(intT i=0; i<n; ++i) {
      cells[i].init();
    }
    numCells = 0;}

  inline void initRand() {
    srand(time(NULL));
    for (intT i = 0; i < dim; i++) {
      if(noRandom) randInt[i] = rands[i] % range + 1;
      else randInt[i] = rand() % range + 1;}
  }

  grid(intT nn, pointT pMinn, pointPairT RR): pMin(pMinn), R(RR), n(nn) {
    initRand();
    cells = newA(cellT, n);
    resetCells();}
  ~grid() {free(cells);}

  void addTable(tableT* tablee) {table = tablee;}
  void updateR(pointPairT RR) {R = RR;}
  pointPairT getR() {return R;}
  floatT r() {return R.dist*(floatT)gridSizeRatio;}

  inline int compareCoordinate(floatT* x1, floatT* x2) {
    for(int i=0; i<dim; ++i) {
      intT xx1 = (intT) floor((x1[i]-pMin[i])/r());
      intT xx2 = (intT) floor((x2[i]-pMin[i])/r());
      if (xx1 != xx2) {
        if (xx1 > xx2) return 1;
        else return -1;
      }}
    return 0;}

  template<class func>
  void cellNeighborhood(cellT* center, func* f) {
    if (center->isEmpty()) {
      cout << "error, finding neighbor for empty cell, abort"<< endl; abort();}
    intT numNeighbor = 3;//3**dim
    for(int i=1; i<dim; ++i) numNeighbor*=3;
    pointT p = pointT();
    cellT bait = cellT(p);
    for(int i=0; i<numNeighbor; ++i) {
      int ii = i;
      for(int j=0; j<dim; ++j) {
        bait.coordP.x[j] = center->P[0][j]+((floatT)(ii%3)-1)*r()*1.00000001;//numerical stability
        ii /= 3;}
      cellT* nbr = table->find(&bait);
      if (!nbr->isEmpty()) {
        if (f->checkComplete(nbr)) break;}
    }
  }

  intT reInsertSerial(iPointT* P, intT nn, bool useTree=false) {
    resetCells();
    auto pLess = [&] (iPointT a, iPointT b) {
                   return pointGridCmp<dim, iPointT>(a, b, pMin, r());
                 };
    quickSortSerial(P, nn, pLess);
    intT c = 0;
    cells[c].P = P;
    for (intT i=1; i<nn; ++i) {
      if (compareCoordinate(P[i].x(), P[i-1].x())) {
        //next cell
        cells[c].numPoints = &P[i]-cells[c].P;
        cells[c].computeCoord(pMin, r());
        c++;
        cells[c].P = &P[i];
      }
    }
    cells[c].numPoints = &P[nn-1]-cells[c].P+1;
    cells[c].computeCoord(pMin, r());
    c++;
    numCells = c;
    if (!useTree) {
      table->setActive(min(n,numCells*3));
      table->clear();
      for (intT i=0; i<numCells; ++i) {
        table->insert(&cells[i]);}
    }
    return numCells;
  }

  intT reInsert(iPointT* P, intT nn, bool useTree=false) {
    resetCells();
    auto pLess = [&] (iPointT a, iPointT b) {
                   return pointGridCmp<dim, iPointT>(a, b, pMin, r());
                 };
    sampleSort(P, nn, pLess);

    int* flag = newA(int, nn+1);
    flag[0] = 1;
    par_for(intT i=1; i<nn; ++i) {
      if (compareCoordinate(P[i].x(), P[i-1].x())) flag[i] = 1;
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
      cells[i].computeCoord(pMin, r());
    }
    cells[numCells-1].numPoints = &P[nn] - cells[numCells-1].P;
    cells[numCells-1].computeCoord(pMin, r());

    if (!useTree) {
      table->setActive(min(n,numCells*3));
      table->clear();
      par_for(intT i=0; i<numCells; ++i) {
        table->insert(&cells[i]);
      }
    }
    free(flag);
    return numCells;
  }

  inline bool findConflict(pointT pp) {
    cellT* target = table->find(new cellT(pp));
    if (target->isEmpty()) {
      cout << "error, findConflict empty target, abort" << endl; abort();}
    auto finder = conflictFinder<dim>(pp, r());
    cellNeighborhood(target, &finder);
    return finder.hasConflict();
  }

  inline pointPairT findClosest(pointT pp) {
    cellT* target = table->find(new cellT(pp));
    if (target->isEmpty()) {
      cout << "error, findClosest empty target, abort" << endl; abort();}
    auto finder = closestFinder<dim>(pp);
    cellNeighborhood(target, &finder);
    return finder.getClosest();}

  inline bool findConflictKdTree(pointT pp, kdTree<dim,cellT>* tree) {
    iPointT p = iPointT(pp);
    cellT bait = cellT(&p,1);
    bait.computeCoord(pMin, r());
    auto finder = conflictFinder<dim>(pp, r());
    tree->rangeNeighbor(&bait, sqrt(2)*r(), &finder);
    return finder.hasConflict();}

  inline pointPairT findClosestKdTree(pointT pp, kdTree<dim, cellT>* tree) {
    iPointT p = iPointT(pp);
    cellT bait = cellT(&p,1);
    bait.computeCoord(pMin, r());
    auto finder = closestFinder<dim>(pp);
    tree->rangeNeighbor(&bait, sqrt(2)*r(), &finder);
    return finder.getClosest();}

  pair<intT, pointPairT> checkInsertSerial(pointT* P, intT s, intT e, bool useTree=false) {
    kdTree<dim, cellT>* tree;
    if (useTree) tree = new kdTree<dim, cellT>(cells, numCells, false);
    for(intT i=s; i<e; ++i) {
      bool hasConflict;
      if (useTree) hasConflict = findConflictKdTree(P[i], tree);
      else hasConflict = findConflict(P[i]);
      if(hasConflict) {
        if (useTree) {
          auto RR = make_pair(i, findClosestKdTree(P[i], tree));
          delete(tree);
          return RR;
        } else {
          return make_pair(i, findClosest(P[i]));
        }
      }
    }
    if (useTree) delete(tree);
    return make_pair(e, pointPairT());
  }

  pair<intT, pointPairT> checkInsert(pointT* P, intT s, intT e, bool useTree=false) {
    kdTree<dim, cellT>* tree;
    if (useTree) tree = new kdTree<dim, cellT>(cells, numCells, true);
    intT bucketSize = min(e-s, (intT)4000);
    intT numBuckets = (e-s)/bucketSize;
    if ((e-s)%bucketSize > 0) numBuckets++;
    intT* bucket = newA(intT, bucketSize+1);

    for (intT b=0; b<numBuckets; ++b) {
      intT ss = s+b*bucketSize;
      intT ee = ss+min(bucketSize, e-ss);
      par_for(intT i=ss; i<ee; ++i) {
        bool hasConflict;
        if (useTree) hasConflict = findConflictKdTree(P[i], tree);
        else hasConflict = findConflict(P[i]);
        if (hasConflict) bucket[i-ss] = 1;
        else bucket[i-ss] = 0;
      }
      intT numBad = sequence::prefixSum(bucket,0,ee-ss);
      bucket[ee-ss] = numBad;

      if (numBad > 0) {
        intT ii = -1;
        par_for(intT i=0; i<ee-ss; ++i) {
          if (bucket[i]==0 && bucket[i]!=bucket[i+1]) ii = i;
        }
        free(bucket);
        if (useTree) {
          auto RR = make_pair(ss+ii, findClosestKdTree(P[ss+ii], tree));
          delete(tree);
          return RR;
        } else return make_pair(ss+ii, findClosest(P[ss+ii]));
      }
    }
    free(bucket);
    if (useTree) delete(tree);
    return make_pair(e, pointPairT());
  }
};

template<int dim>
struct cellHash {
  typedef grid<dim> gridT;
  typedef cell<dim>* eType;
  typedef cell<dim>* kType;
  gridT* grids;
  cellHash(gridT* gridss): grids(gridss) {};
  eType empty() {return new cell<dim>();}
  kType getKey(eType v) {return v;}
  uintT hash(kType c) {
    if (c->isEmpty()) {// todo, remove
      cout << "error, hashing empty cell, abort" << endl; abort();}
    return grids->hash(c->coordinate());}
  int cmp(kType c1, kType c2) {
    if (c1->isEmpty() || c2->isEmpty()) return 1;
    return grids->compareCoordinate(c1->coordinate(), c2->coordinate());}
  bool replaceQ(eType c1, eType c2) {return 0;}
};

// *************************************************************
//    Incremental algorithm
// *************************************************************

/**
 * Computes the closest pair of P using the randomized incremental algorithm.
 * @param P a point array.
 * @param n length of PIn.
 * @param serial whether to run serially.
 * @return the closest pair
 */
template<int dim>
pointPair<dim> incremental(point<dim>* P, intT n, bool serial=false) {
  typedef point<dim> pointT;
  typedef iPoint<dim> iPointT;
  typedef pointPair<dim> pointPairT;
  static const bool verbose = true;
  static const bool noRandom = false;
  bool useTree = false;
  if(dim>=5) useTree = true;

  timing t0; t0.start();
  double tInsert = 0;
  double tCheck = 0;

  static intT cutoff = 1000;
  pointPairT R;
  if(serial) R = bruteForceSerial<dim>(P, min(n, cutoff));
  else R = bruteForceParallel<dim>(P, min(n, cutoff));
  if (n <= cutoff) return R;

  if(!noRandom) {
    //permutation
    if(verbose) cout << "permuting points" << endl;
    randPerm(P, n);
  }

  pointT pMin;
  if(serial) pMin = pMinSerial(P, n);
  else pMin = pMinParallel(P, n);

  auto grids = grid<dim>(n, pMin, R);
  auto hash = cellHash<dim>(&grids);
  Table<cellHash<dim>,intT>* table;
  if (useTree) table = new Table<cellHash<dim>,intT>(1, hash);
  else table = new Table<cellHash<dim>,intT>(n, hash, 4);
  grids.addTable(table);

  iPointT *PP = newA(iPointT, n);
  par_for(intT i=0; i<n; ++i) {
    PP[i]=iPointT(P[i], i);
  }

  grids.updateR(R);
  intT nn = 0;
  intT rounds = 0;
  while (nn < n) {
    rounds++;
    t0.start();
    intT toInsert;
    if (nn != 0) {
      toInsert = min(n-nn, nn);
    } else {
      toInsert = min(n-1000, (intT)1000);
    }

    if(serial) grids.reInsertSerial(PP, nn+toInsert, useTree);
    else grids.reInsert(PP, nn+toInsert, useTree);
    tInsert += t0.next();

    pair<intT, pointPairT> insert;
    if(serial) insert = grids.checkInsertSerial(P, nn, nn+toInsert, useTree);
    else insert = grids.checkInsert(P, nn, nn+toInsert, useTree);

    intT inserted = insert.first; pointPairT newR = insert.second;
    if (!newR.isEmpty()) {
      if(verbose) cout << "r = " << newR.dist;
      grids.updateR(newR);
      auto pLess = [&] (iPointT a, iPointT b) {
                     return (a.idx() == a.idx()) ? false :
                       (a.idx() < b.idx() ? true : false);};

      if(serial) quickSortSerial(P, nn+toInsert, pLess);
      else sampleSort(P, nn+toInsert, pLess);
    } else if (newR.dist == 0) {
      grids.updateR(newR);
      t0.stop();
      break;
    } else {
      if(verbose)cout << "r = " << grids.r();
    }
    //if(verbose) cout << "round inserted " << inserted-nn << "/" << toInsert << endl;
    nn = inserted;
    if(verbose) cout << ", total inserted " << nn << "/" << n << endl;
    tCheck += t0.stop();
  }
  cout << "insert-time = " << tInsert << endl;
  cout << "check-time = " << tCheck << endl;
  cout << rounds << " rounds" << endl;
  free(PP);
  delete(table);
  return grids.getR();
}

// *************************************************************
//    DRIVER
// *************************************************************

/**
 * Computes the closest pair of P using the randomized incremental algorithm.
 * @param P a point array.
 * @param n length of PIn.
 * @return the closest pair
 */
template<int dim>
pair<point<dim>, point<dim>> closestPair(point<dim>* P, intT n) {
  static const bool serial = false;
  cout << "closestPair of " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  pointPair<dim> r = incremental<dim>(P, n, serial);
  cout << "incremental " << r.u << ", " << r.v << ", dist " << r.dist << endl << endl;

  return make_pair(r.u, r.v);
}

template pair<point<2>, point<2>> closestPair<2>(point<2>*, intT);
template pair<point<3>, point<3>> closestPair<3>(point<3>*, intT);
template pair<point<4>, point<4>> closestPair<4>(point<4>*, intT);
template pair<point<5>, point<5>> closestPair<5>(point<5>*, intT);
template pair<point<6>, point<6>> closestPair<6>(point<6>*, intT);
template pair<point<7>, point<7>> closestPair<7>(point<7>*, intT);
template pair<point<8>, point<8>> closestPair<8>(point<8>*, intT);
template pair<point<9>, point<9>> closestPair<9>(point<9>*, intT);
