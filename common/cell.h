// Copyright (c) 2020 Yiqiu Wang
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

#ifndef CELL_H
#define CELL_H

#include "grid.h"
#include "pbbs/ndHash.h"

#ifdef USEJEMALLOC
#include<jemalloc/jemalloc.h>
#define jeNewA(__E,__n) (__E*) je_custom_prefix_malloc((__n)*sizeof(__E))
#define jeFree(__E) je_custom_prefix_free(__E)
#endif

template<int dim, class objT> struct grid;
template<int dim, class objT> struct cellHash;

/**
 *  A cell class that represents each box in the grid.
 */
template<int dim, class objT>
struct cell {
  typedef double floatT;
  typedef point<dim> geoPointT;
  typedef cell<dim, objT> cellT;
  typedef Table<cellHash<dim, objT>,intT> tableT;
  typedef grid<dim, objT> gridT;

  static const intT defaultCapacity = 10;
  static const intT resizeFactor = 2;

  objT *P=NULL;
  geoPointT coordP;
  intT numPoints;
  intT numErased;
  intT capacity;

  /**
  *   Constructor 1, that initializes an empty cell.
  */
  cell(): numPoints(0), numErased(0), capacity(0) {};

  /**
  *   Constructor 2, that initializes cell with external point array.
  *   @param PP external point array.
  *   @param nn size of PP.
  */
  cell(objT* PP, intT nn): P(PP), numPoints(nn), numErased(0), capacity(0){};

  /**
  *   Constructor 3, that only initializes the coordinate of the cell base on an external point.
  *   @param coordPP external point.
  */
  cell(geoPointT coordPP): coordP(coordPP), numPoints(0), numErased(0), capacity(0){}//only used as bait
  cell(floatT* coordIn): numPoints(0), numErased(0), capacity(0){
    for(intT i=0; i<dim; ++i) coordP.updateCoordinate(i, coordIn[i]);
  }//only used as bait

  inline void init() {
    capacity = 0;
    numPoints = 0;
    numErased = 0;}

  ~cell() {
    del();
  }

  inline void del() {
    if(P&&capacity>0) {
#ifdef USEJEMALLOC
      jeFree(P);
#else
      free(P);
#endif
    }
  }

  /**
  *   Resizes the cell to the size wanted. //todo see if x2 improves perf
  *   @param want size wanted.
  */
  void resize(intT want) {
    auto oldCapacity = capacity;
    capacity = max(want, numPoints);
    if (capacity <= 0) capacity = defaultCapacity;
#ifdef USEJEMALLOC
    auto newP = jeNewA(objT, capacity);
#else
    auto newP = newA(objT, capacity);
#endif
    if (numPoints > 2000) {
      par_for(intT i=0; i<numPoints; ++i) {newP[i]=P[i];}
    }
    else {for(intT i=0; i<numPoints; ++i) newP[i] = P[i];}

    if (capacity-numPoints > 2000) {
      par_for(intT i=numPoints; i<capacity; ++i) {newP[i].setEmpty();}
    }
    else {for(intT i=numPoints; i<capacity; ++i) newP[i].setEmpty();}
    swap(newP, P);
    if (oldCapacity>0) {
#ifdef USEJEMALLOC
      jeFree(newP);
#else
      free(newP);
#endif
    }
  }

  /**
  *   Compacts cell to get rid of deleted entries if real size is < 0.25 allocated.
  *   Halfs the memory usage.
  */
  inline void compact() {
    if (actualSize() <= 0) {
      if(capacity >0) {
#ifdef USEJEMALLOC
        jeFree(P);
#else
        free(P);
#endif
      }
      init();
      return;}

    if (size() < actualSize()*4) return;

    auto oldCapacity = capacity;
    capacity = max(capacity/2, (intT)2);
#ifdef USEJEMALLOC
    auto newP = jeNewA(objT, capacity);
#else
    auto newP = newA(objT, capacity);
#endif
    intT ii=0;
    for(intT i=0; i<numPoints; ++i) {//todo consider make parallel
      if(!P[i].isEmpty()) newP[ii++] = P[i];}
    numPoints = ii;
    for(intT i=numPoints; i<capacity; ++i) {
      newP[i].setEmpty();}
    numErased = 0;
    swap(newP, P);
    if (oldCapacity>0) {
#ifdef USEJEMALLOC
      jeFree(newP);
#else
      free(newP);
#endif
    }
  }

  /**
  *   Determine if cell contains a duplicate entry.
  *   @param p the entry whose duplicate is to be detected.
  *   @return yes or no.
  */
  inline bool dupDetector(objT p) {
    for (intT i=0; i<size(); ++i) {
      if (!P[i].isEmpty() && samePoint<dim>(P[i].coordinate(), p.coordinate())) {
        return true;}
    }
    return false;}

  /**
  *   Insert p if no duplicate of p is in the cell.
  *   @param p the point to be inserted.
  *   @return whether insertion happened.
  */
  bool insert(objT p) {
    if (dupDetector(p)) return false;
    if(numPoints+1>=capacity || capacity <= 0) resize((numPoints+1)*resizeFactor);
    P[numPoints++] = p;
    return true;}

  inline void incErased() {numErased++;}

  /**
  *   Delete p from the cell.
  *   @param p the point to be inserted.
  *   @return whether deletion happened.
  */
  bool erase(objT p) {
    for (intT i=0; i<size(); ++i) {
      if (!P[i].isEmpty() && samePoint<dim>(P[i].coordinate(), p.coordinate())) {
        P[i].setEmpty();
        incErased();
        compact();
        return true;}
    }
    return false;
  }

  /**
  *   Computes the coordinate of the cell base on P, and assign to coordP.
  *   The coordinate is suppose to be the center of the cell.
  *   @param pMin the global point minimum.
  *   @param r the grid size.
  */
  void computeCoord(geoPointT pMin, double r) {
    for(int i=0; i<dim; ++i) {
      coordP.x[i] = r/2+pMin[i]+(double)floor((P[0][i]-pMin[i])/r)*r;}
  }

  /**
  *   Computes the coordinate of the cell base on coordP, and update coordP.
  *   The coordinate is suppose to be the center of the cell.
  *   @param pMin the global point minimum.
  *   @param r the grid size.
  */
  void coordCenter(geoPointT pMin, double r) {
    for(int i=0; i<dim; ++i) {
      coordP.x[i] = r/2+pMin[i]+(double)floor((coordP[i]-pMin[i])/r)*r;}
  }

  /**
  *   The number of inserted points, not considering deletion.
  *   @return number.
  */
  inline intT size() {return numPoints;}

  /**
  *   The number of inserted points (&& not deleted yet).
  *   @return number.
  */
  inline intT actualSize() {return numPoints-numErased;}

  /**
  *   Whether the cell is valid (initialized). It has nothing to do with the number of inserted points.
  *   @return yes or no.
  */
  inline bool isEmpty() {return coordP.isEmpty();}

  floatT *coordinate() {
    if (isEmpty()) return NULL;
    else return coordP.x;
  }

  floatT coordinate(int i) {return coordP.x[i];}

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

  void printPoints() {
    for (intT i=0; i<size(); ++i) {
      cout << P[i];
      if (P[i].isEmpty()) cout << "x";
      cout << " ";
    }
    cout << endl;
  }

  template<class func>
  void pointMap(func& f) {
    for(intT i=0; i<numPoints; ++i) {
      if(!P[i].isEmpty()) f(P[i]);}
  }
};

/**
  *   Hash function for cell, for pbbs/ndHash.h
  */
template<int dim, class objT>
struct cellHash {
  typedef double floatT;
  typedef cell<dim, objT> cellT;
  typedef hashFloatToCell<dim> hashFunc;
  typedef cellT* eType;
  typedef cellT* kType;

  hashFunc* hashF;
  cellT* e;

  cellHash(hashFunc* hashFF):hashF(hashFF) {
    e = new cellT();}

  ~cellHash() {}

  eType empty() {return e;}

  kType getKey(eType v) {return v;}

  uintT hash(kType c) {
    return hashF->hash(c->coordinate());
  }

  int cmp(kType c1, kType c2) {
    if (c1->isEmpty() || c2->isEmpty()) return 1;
    return hashF->compareCell(c1->coordinate(), c2->coordinate());
  }

  inline int diffCell(floatT* c1, floatT* c2) {return hashF->compareCell(c1, c2);}

  bool replaceQ(eType c1, eType c2) {return 0;}
};

#endif
