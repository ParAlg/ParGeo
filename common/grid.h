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

#ifndef GRID_H
#define GRID_H

#include "cell.h"
#include "geometry.h"
#include "shared.h"
#include "dynamicKdTree.h"
#include "pbbs/sequence.h"
#include "pbbs/ndHash.h"
#include "pbbs/sampleSort.h"
#include "pbbs/quickSort.h"
#include "pbbs/parallel.h"

#ifdef USEJEMALLOC
#include<jemalloc/jemalloc.h>
#define jeNewA(__E,__n) (__E*) je_custom_prefix_malloc((__n)*sizeof(__E))
#define jeFree(__E) je_custom_prefix_free(__E)
#endif

//a less comparator based on grid
template<int dim, class pointT, class geoPointT>
inline bool pointGridCmp(pointT p1, pointT p2, geoPointT pMin, floatT r) {
  for(int i=0; i<dim; ++i) {
    intT xx1 = (intT) floor((p1[i]-pMin[i])/r);
    intT xx2 = (intT) floor((p2[i]-pMin[i])/r);
    if (xx1 != xx2) {
      if (xx1 > xx2) return false;
      else return true;}
  }
  return false;
}

/**
  *   A grid class, that puts dim-dimensional axis-aligned box cells on a point set.
  */
template<int dim, class objT>
struct grid {
  typedef grid<dim, objT> gridT;
  typedef double floatT;
  typedef point<dim> geoPointT;
  typedef cell<dim, objT> cellT;
  typedef hashFloatToCell<dim> cellHashT;
  typedef Table<cellHash<dim, objT>,intT> tableT;
  typedef Table<aFloatHash<dim, objT>,intT> objTableT;
  typedef dynamicKdTree<dim, cellT> treeT;
  typedef pointPair<dim> pointPairT;

  static const bool noRandom = true;

  bool useTree = false;

  //pivotT v;
  floatT r;
  geoPointT pMin;
  cellT* cells;
  intT numCells, cellCapacity;
  cellHashT* myHash=NULL;// generic hash function
  tableT* table=NULL;
  treeT* tree=NULL;
  intT totalPoints;

  inline void resetCells() {
    parallel_for(0, cellCapacity, [&](intT i) {cells[i].init();});
    numCells = 0;}

  /**
  *   Grid constructor.
  *   @param cellMax projected maximum number of points inserted.
  *   @param pMinn global coordinate minimum.
  *   @param r box cell size.
  */
  grid(intT cellMax, geoPointT pMinn, floatT rr):
    pMin(pMinn), cellCapacity(cellMax), totalPoints(0), r(rr) {
    cells = newA(cellT, cellCapacity);
    resetCells();
    myHash = new cellHashT(pMinn, r);
    table = new tableT(cellMax*2, cellHash<dim, objT>(myHash));//todo load

    if(dim>=5) useTree = true;
    if(useTree) tree = new treeT();
  }

  ~grid() {
    parallel_for(0, numCells, [&](intT i) {cells[i].del();});
    free(cells);
    if(myHash) delete myHash;
    if(table) {
      table->del();
      delete table;}
    if(useTree && tree) delete tree;
  }

  /**
  *   Find a cell base on an arbitrary coordinate that falls in the cell.
  *   @param coord coordinate of length dim.
  *   @return reference to the cell, if not found then return empty cell.
  */
  inline cellT* getCell(floatT* coord) {
    cellT bait = cellT(geoPointT(coord));
    cellT* found = table->find(&bait);
    return found;}

  //number of points in all cells
  inline intT size() {
    return totalPoints;
  }

  /**
  *   Serially, computes upper bound of the sum of sizes of the neighborhoods of a set of points.
  *   @param PP point array.
  *   @param nn length of PP.
  *   @return total size.
  */
  intT nghSizeSerial(objT* PP, intT nn) {
    intT total = 0;
    auto fSize = [&](cellT& c) {
                   total += c.actualSize();
                   return false;};
    for(intT i=0; i<nn; ++i) {
      nghCellMap(PP[i].coordinate(), fSize);}
    return total;}

  /**
  *   In parallel, computes upper bound of the sum of sizes of the neighborhoods of a set of points.
  *   @param ...
  *   @param flag auxiliary memory of length nn+1.
  *   @return ...
  */
  intT nghSizeParallel(objT* PP, intT nn, intT* flag=NULL) {
    bool freeFlag=false;
    if(!flag) {
      flag = newA(intT, nn+1);
      freeFlag=true;}
    parallel_for(0, nn, [&](intT i) {
	flag[i] = 0;
	auto fSize = [&](cellT& c) {
	  flag[i] += c.actualSize();
	  return false;};
	nghCellMap(PP[i].coordinate(), fSize);
      });
    intT total = sequence::prefixSum(flag, 0, nn);
    if(freeFlag) free(flag);
    return total;}

  /**
  *   Map across the neighborhood of a cell that contains a given coordinate (include self), and perform actions on those cells.
  *   The function automatically switch to kdTree for higher dimensions.
  *   @param center the coordinate that the cell contains.
  *   @param f a lambda function that takes in a point.
  */
  template<class func>
  inline void nghPointMap(floatT* center, func& f) {
    if(useTree) {
      geoPointT p = geoPointT();
      cellT bait = cellT(p);
      for(int j=0; j<dim; ++j) bait.coordP.x[j] = center[j];
      bait.coordCenter(pMin, r);
      auto fStop = [&](){return false;};
      auto fWrap = [&](cellT* nbr) {
                     if (!nbr->isEmpty()
                         && nbr->actualSize()>0
                         && nbr->nextTo(bait, r)
                         ) {
                       for(intT jj=0;jj<nbr->size();++jj) {
                         if(!nbr->P[jj].isEmpty()) {
                           if(f(nbr->P[jj])) return true;}//stop iteration
                       }
                     }
                     return false;};//todo, optimize
      tree->rangeNeighbor(&bait, r*sqrt(dim), fWrap, fStop);
    } else {
      intT numNeighbor = 3;//3**dim
      for(int i=1; i<dim; ++i) numNeighbor*=3;
      geoPointT p = geoPointT();
      cellT bait = cellT(p);
      for(int i=0; i<numNeighbor; ++i) {
        int ii = i;
        for(int j=0; j<dim; ++j) {
          bait.coordP.x[j] = center[j]+((floatT)(ii%3)-1)*r*1.000001;
          ii /= 3;}
        cellT* nbr = table->find(&bait);
        if (!nbr->isEmpty() && nbr->actualSize()>0) {
          for(intT jj=0;jj<nbr->size();++jj) {
            if(!nbr->P[jj].isEmpty()) {
              if(f(nbr->P[jj])) break;}
          }
        }}
    }
  }

  /**
  *   Map across the neighborhood cells of a cell that contains a given coordinate (include self), and perform actions on those cells.
  *   The function automatically switch to kdTree for higher dimensions.
  *   @param center the coordinate that the cell contains.
  *   @param f a lambda function that takes in reference of a cell.
  */
  template<class func>
  inline void nghCellMap(floatT* center, func& f) {
    if(useTree) {
      geoPointT p = geoPointT();
      cellT bait = cellT(p);
      for(int j=0; j<dim; ++j) bait.coordP.x[j] = center[j];
      bait.coordCenter(pMin, r);
      auto fStop = [&](){return false;};
      auto fWrap = [&](cellT* cell){
                     if(!cell->isEmpty()
                        && cell->nextTo(bait, r))
                       return f(*cell);
                     return false;
                   };
      tree->rangeNeighbor(&bait, r*sqrt(dim), fWrap, fStop);
    } else {
      intT numNeighbor = 3;//3**dim
      for(int i=1; i<dim; ++i) numNeighbor*=3;
      geoPointT p = geoPointT();
      cellT bait = cellT(p);
      for(int i=0; i<numNeighbor; ++i) {
        int ii = i;
        for(int j=0; j<dim; ++j) {
          bait.coordP.x[j] = center[j]+((floatT)(ii%3)-1)*r*1.000001;
          ii /= 3;}
        cellT* nbr = table->find(&bait);
        if (!nbr->isEmpty()) {
          if(f(*nbr)) break;}
      }
    }
  }

  /**
  *   Serially, map across set of points pNgh in the neighborhood of an input set of points. Apply a function on a specified subset of pNgh.
  *   The function automatically switch to kdTree for higher dimensions.
  *   @param PP input array of points.
  *   @param nn length of PP.
  *   @param filter a lambda function that takes in a point of pNgh, returns if want to process the point.
  *   @param map a lambda function that takes in an index and a point of pNgh, and processes the point. The index is out of the the processed points of pNgh only.
  *   @param flag auxiliary array of size nn+1.
  *   @return the size of the subset of pNgh processed.
  */
  template<class filterFuncT, class mapFuncT>
  intT nghPointMapSerial(objT* PP, intT nn, filterFuncT filter, mapFuncT map) {
    intT nNbrs = 0;
    auto fSum = [&](cellT& c) {nNbrs+=c.actualSize();return false;};
    for(intT i=0; i<nn; ++i) {
      nghCellMap(PP[i].coordinate(), fSum);
    }

    if(nNbrs<=0) return 0;

    objT** pointers = newA(objT*, nNbrs);
    intT ii=0;
    auto fGetPt = [&](objT& q) {
                    if(filter(q)) pointers[ii++]=&q;
                    return false;};
    for(intT i=0; i<nn; ++i) {
      nghPointMap(PP[i].coordinate(), fGetPt);}

    if(ii<=0) {
      free(pointers);
      return 0;
    }

    if (noRandom) {
      auto pLess = [&](objT *a, objT *b) {return (*a)<(*b);};
      quickSortSerial(pointers, ii, pLess);
    } else {
      auto ppLess = [&](objT *a, objT *b) {return a<b;};
      quickSortSerial(pointers, ii, ppLess);
    }

    intT iii=0;
    for(intT i=0; i<ii; ++i) {
      if(i>0 && pointers[i]==pointers[i-1]) continue;
      map(iii++, *(pointers[i]));
    }
    free(pointers);
    return iii;
  }

  /**
  *   In parallel, map across set of points pNgh in the neighborhood of an input set of points. Apply a function on a specified subset of pNgh.
  *   The function automatically switch to kdTree for higher dimensions.
  *   @param ...
  *   @return ...
  */
  template<class filterFuncT, class mapFuncT>
  intT nghPointMapParallel(objT* PP, intT nn, filterFuncT filter, mapFuncT map, intT* flag=NULL) {
    bool freeFlag=false;
    if(!flag) {
      flag=newA(intT, nn+1);
      freeFlag=true;}

    parallel_for(0, nn, [&](intT i) {
	flag[i] = 0;
	auto fSum = [&](cellT& c) {flag[i]+=c.actualSize();return false;};
	nghCellMap(PP[i].coordinate(), fSum);
      });
    intT nNbrs = sequence::prefixSum(flag, 0, nn);
    if(nNbrs <= 0) {
      if(freeFlag) {
        free(flag);}
      return 0;
    }
    flag[nn] = nNbrs;

    objT** pointers = newA(objT*, nNbrs);
    parallel_for(0, nn, [&](intT i) {
	intT ii=0;
	auto fGetPt = [&](objT& q) {
	  if(filter(q)) pointers[flag[i]+(ii++)]=&q;
	  else pointers[flag[i]+(ii++)]=NULL;
	  return false;};
	nghPointMap(PP[i].coordinate(), fGetPt);
      });

    if (noRandom) {
      auto pLess = [&](objT *a, objT *b) {
                      if(!a) return false;//put NULL at end
                      if(!b) return true;//put NULL at end
                      return (*a)<(*b);};
      sampleSort(pointers, nNbrs, pLess);
    } else {
      auto ppLess = [&](objT *a, objT *b) {
                      if(!a) return false;//put NULL at end
                      if(!b) return true;//put NULL at end
                      return a<b;};
      sampleSort(pointers, nNbrs, ppLess);
    }

    intT* flag2=newA(intT, nNbrs+1);

    if(pointers[0]) flag2[0] = 1;
    else flag2[0] = 0;
    parallel_for(1, nNbrs, [&](intT i) {
	if(pointers[i] && pointers[i] != pointers[i-1]) flag2[i] = 1;
	else flag2[i] = 0;
      });
    intT nPts = sequence::prefixSum(flag2, 0, nNbrs);
    flag2[nNbrs] = nPts;

    if(nPts <= 0) {
      free(pointers);
      if(freeFlag) free(flag);
      free(flag2);
      return 0;
    }

    parallel_for(0, nNbrs, [&](intT i) {
	if(flag2[i] != flag2[i+1]) map(flag2[i], *(pointers[i]));
      });
    free(pointers);
    if(freeFlag) free(flag);
    free(flag2);
    return nPts;
  }

  /**
  *   Serial, hash-based nghPointMap.
  *   @param ...
  *   @param sorted if PP is sorted based on cell.
  *   @return ...
  */
  template<class filterFuncT, class mapFuncT>
  intT nghPointMapHashSerial(objT* PP, intT nn, filterFuncT filter, mapFuncT map, bool sorted=false) {
    if(!sorted) {
      auto pLess = [&] (objT a, objT b) {
                     return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
      quickSortSerial(PP, nn, pLess);}

    auto pTable = objTableT(nghSizeSerial(PP, nn), aFloatHash<dim, objT>(myHash));
    auto fInsert = [&] (objT& p) {
                     if(filter(p)) pTable.insert(&p);
                     return false;};
    for(intT i=0; i<nn; ++i) {
      if (i==0 || table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {//short circuit
        nghPointMap(PP[i].coordinate(), fInsert);}
    }

    auto P0 = pTable.entries();
    pTable.del();
    for(intT i=0; i<P0.n; ++i) {
      map(i, *P0.A[i]);}
    free(P0.A);
    return P0.n;
  }

  /**
  *   Parallel, hash-based nghPointMap.
  *   @param ...
  *   @param sorted if PP is sorted based on cell.
  *   @return ...
  */
  template<class filterFuncT, class mapFuncT>
  intT nghPointMapHashParallel(objT* PP, intT nn, filterFuncT filter, mapFuncT map, bool sorted=false, intT* flag=NULL) {
    if(!sorted) {
      auto pLess = [&] (objT a, objT b) {
                     return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
      sampleSort(PP, nn, pLess);}

    auto pTable = objTableT(nghSizeParallel(PP, nn, flag), aFloatHash<dim, objT>(myHash));
    auto fInsert = [&] (objT& p) {
                     if(filter(p)) pTable.insert(&p);
                     return false;};
    parallel_for(0, nn, [&](intT i) {
	if (i==0 || table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {//short circuit
	  nghPointMap(PP[i].coordinate(), fInsert);}
      });

    auto P0 = pTable.entries();
    pTable.del();
    parallel_for(0, P0.n, [&](intT i) {
	map(i, *P0.A[i]);});
    free(P0.A);
    return P0.n;
  }

  template<class mapFuncT>
  inline void allCellMap(mapFuncT f) {
    parallel_for(0, numCells, [&](intT i) {
	if (!cells[i].isEmpty()) f(&cells[i]);});
  }

  /**
  *   Serially, map across all points in the grid.
  *   @param f map function that takes in index and point.
  *   @return total number of points
  */
  template<class mapFuncT>
  inline intT allPointMapSerial(mapFuncT f) {
    intT count = 0;
    for(intT i=0; i<numCells; ++i) {
      auto c = &cells[i];
      for (intT j=0; j<c->size(); ++j) {
        if (!c->P[j].isEmpty()) {
          f(count++, c->P[j]);
        }}
    }
    return count;
  }

  /**
  *   In parallel, map across all points in the grid.
  *   @param ...
  *   @param counts auxiliary memory of size numCells+1
  *   @return ...
  */
  template<class mapFuncT>
  inline intT allPointMapParallel(mapFuncT f, intT* counts=NULL) {
    bool freeCounts=false;
    if (!counts) {
      counts = newA(intT, numCells+1);
      freeCounts = true;}
    parallel_for(0, numCells, [&](intT i) {
	counts[i]=cells[i].actualSize();});
    intT nPoints = sequence::prefixSum(counts, 0, numCells);
    parallel_for(0, numCells, [&](intT i) {
	auto c = &cells[i];
	intT jj=0;
	for (intT j=0; j<c->size(); ++j) {
	  if (!c->P[j].isEmpty()) {
	    f(counts[i]+jj, c->P[j]);
	    jj++;}
	}});
    if(freeCounts) free(counts);
    return nPoints;
  }

  //true if point nghood is empty
  inline bool isSparse(objT pp) {
    intT numPoints=0;
    auto fSparse = [&] (cellT& c) {
                     numPoints += c.actualSize();return false;};
    nghCellMap(pp.coordinate(), fSparse);
    return numPoints <=1;
  }

  /**
  *   Serially, erase a set of points from the grid, and optionally return their de-duplicated neighbors.
  *   @param PP input point array.
  *   @param nn size of PP.
  *   @param returnNeighbors whether to return the de-duplicated neighbors of PP.
  *   @param P1 auxiliary memory of size nn.
  *   @return sequence of de-duplicated neighbors.
  */
  _seq<objT> eraseSerial(objT* PP, intT nn, bool returnNeighbors=false, objT* P1=NULL) {
    if (nn<=0) return _seq<objT>();
    auto pLess = [&] (objT a, objT b) {
                   return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
    quickSortSerial(PP, nn, pLess);
    intT ii=0;
    cellT* c;
    auto bait = cellT(PP[0].coordinate());
    c = table->find(&bait);
    if(!c->isEmpty()) {
      if(c->erase(PP[0])) totalPoints--;
    }
    for (intT i=1; i<nn; ++i) {
      if (table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {
        auto bait = cellT(PP[i].coordinate());
        c = table->find(&bait);}
      if(!c->isEmpty()) {
        if(c->erase(PP[i])) totalPoints--;
      }
    }

    if (returnNeighbors) {
      if(!P1) P1 = newA(objT, nghSizeSerial(PP, nn));
      auto fGet = [&](intT i, objT p) {P1[i] = p;};
      intT total = nghPointMapSerial(PP, nn, [&](objT p){return true;}, fGet);
      //intT total = nghPointMapHashSerial(PP, nn, [&](pointT p){return true;}, fGet, true);
      return _seq<objT>(P1, total);
    }
    return _seq<objT>();
  }

  /**
  *   In parallel, erase a set of points from the grid, and optionally return their de-duplicated neighbors.
  *   @param ...
  *   @param flag auxiliary memory of size 2*(nn+2).//todo check size
  *   @param P1 auxiliary memory of size nn.
  *   @return ...
  */
  _seq<objT> eraseParallel(objT* PP, intT nn, intT* flag=NULL, bool returnNeighbors=false, objT* P1=NULL) {
    if (nn<=0) return _seq<objT>();
    if (nn<=400) return eraseSerial(PP, nn, returnNeighbors, P1);
    bool freeFlag = false;
    if(!flag) {
      flag=newA(intT, 2*(nn+2));
      freeFlag=true;}
    auto pLess = [&] (objT a, objT b) {
                   return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
    sampleSort(PP, nn, pLess);
    parallel_for(0, nn, [&](intT i) {
	if (i==0 || table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {//short circuit
	  auto bait = cellT(PP[i].coordinate());
	  if (!table->find(&bait)->isEmpty())
	    flag[i]=1;
	  else
	    flag[i]=0;
	} else flag[i]=0;
      });
    intT numDelCells = sequence::prefixSum(flag, 0, nn);
    flag[nn] = numDelCells;

    intT* flag2 = flag+nn+1;
    parallel_for(0, nn, [&](intT i) {
	if (flag[i] != flag[i+1]) {
	  auto bait = cellT(PP[i].coordinate());
	  cellT* c = table->find(&bait);
	  auto tmp = c->actualSize();
	  intT j=0;
	  do {
	    c->erase(PP[i+j]);
	    j++;
	  } while (i+j<nn && flag[i+j]==flag[i+j+1]);
	  flag2[i] = tmp - c->actualSize();
	} else flag2[i] = 0;
      });
    totalPoints -= sequence::prefixSum(flag2, 0, nn);

    if (returnNeighbors) {
      if(!P1) P1 = newA(objT, nghSizeParallel(PP, nn, flag));
      auto fGet = [&](intT i, objT p) {P1[i] = p;};
      //intT total = nghPointMapHashParallel(PP, nn, [&](objT p){return true;}, fGet, true, flag);
      intT total = nghPointMapParallel(PP, nn, [&](objT p){return true;}, fGet, flag);
      if(freeFlag) free(flag);
      return _seq<objT>(P1, total);
    }

    if(freeFlag) free(flag);
    return _seq<objT>();
  }

  /**
  *   Serially, insert a set of points to the grid.
  *   @param PP input point array.
  *   @param nn size of PP.
  */
  void insertSerial(objT* PP, intT nn) {
    if (nn==0) return;
    auto pLess = [&] (objT a, objT b) {
                   return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
    quickSortSerial(PP, nn, pLess);
    cellT* c;
    bool oldCell = true;
    auto bait = cellT(PP[0].coordinate());
    c = table->find(&bait);
    if (c->isEmpty()) {
      oldCell = false;
      c = &cells[numCells++];}
    if(c->insert(PP[0])) totalPoints++;
    for (intT i=1; i<nn; ++i) {
      if (table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {
        if (!oldCell) {
          c->computeCoord(pMin, r);
          table->insert(c);
          if(useTree) tree->insert(c, 1);
        }
        oldCell = true;
        auto bait = cellT(PP[i].coordinate());
        c = table->find(&bait);
        if (c->isEmpty()) {
          oldCell = false;
          c = &cells[numCells++];}
      }
      if(c->insert(PP[i])) totalPoints++;;
    }
    if (!oldCell) {
      c->computeCoord(pMin, r);
      table->insert(c);
      if(useTree) tree->insert(c, 1);
    }
    if (numCells > cellCapacity) {
      cout << "error, grid insert exceeded cell capacity, abort()" << endl;abort();}
  }

  /**
  *   In parallel, insert a set of points to the grid.
  *   @param ...
  *   @param flag auxiliary array of size nn*2+2
  *   @param flag auxiliary array of size nn
  */
  void insertParallel(objT* PP, intT nn, intT* flag=NULL, intT* flag3=NULL) {
    if (nn==0) return;
    if (nn<400) insertSerial(PP, nn);
    bool freeFlag=false;
    if (!flag) {
      flag=newA(intT, nn*2+2);
      freeFlag=true;}
    intT* flag2 = flag+nn+1;
    auto pLess = [&] (objT a, objT b) {
                   return pointGridCmp<dim, objT, geoPointT>(a, b, pMin, r);};
    sampleSort(PP, nn, pLess);

    auto bait = cellT(PP[0].coordinate());
    if (!table->find(&bait)->isEmpty()) {
      flag[0] = 1;
      flag2[0] = 0;
    } else {
      flag[0] = 0;
      flag2[0] = 1;}
    parallel_for(1, nn, [&](intT i) {
	if (table->hashStruct.diffCell(PP[i].coordinate(), PP[i-1].coordinate())) {
	  auto bait = cellT(PP[i].coordinate());
	  if (!table->find(&bait)->isEmpty()) {
	    flag[i] = 1;
	    flag2[i] = 0;
	  } else {
	    flag[i] = 0;
	    flag2[i] = 1;}
	} else {
	  flag[i] = 0;
	  flag2[i] = 0;}
      });
    intT numOldCells = sequence::prefixSum(flag, 0, nn);
    flag[nn] = numOldCells;
    intT numNewCells = sequence::prefixSum(flag2, 0, nn);
    flag2[nn] = numNewCells;
    //cout << "num old cells = " << numOldCells << endl;
    //cout << "num new cells = " << numNewCells << endl;

    bool freeFlag3 = false;
    if(!flag3) {
      flag3 = newA(intT, nn);
      freeFlag3 = true;}
    parallel_for(0, nn, [&](intT i) {
	//old cell
	if (flag[i] != flag[i+1]) {
	  auto bait = cellT(PP[i].coordinate());
	  auto c = table->find(&bait );
	  intT tmp = c->actualSize();
	  intT jj=0;
	  do {
	    c->insert(PP[i+jj]);
	    jj++;
	  } while(i+jj<nn
		  && flag[i+jj]==flag[i+jj+1]
		  && flag2[i+jj]==flag2[i+jj+1]);
	  flag3[i] = c->actualSize()-tmp;
	}

	//new cell
	else if (flag2[i] != flag2[i+1]) {
	  auto c = &cells[numCells+flag2[i]];
	  intT tmp = c->actualSize();
	  intT jj=0;
	  do {
	    c->insert(PP[i+jj]);
	    jj++;
	  } while(i+jj<nn
		  && flag[i+jj]==flag[i+jj+1]
		  && flag2[i+jj]==flag2[i+jj+1]);
	  c->computeCoord(pMin, r);
	  flag3[i] = c->actualSize()-tmp;
	  table->insert(c);
	}

	else flag3[i] = 0;
      });

    if(useTree) {
      tree->insert(&cells[numCells], numNewCells);
    }

    numCells += numNewCells;
    totalPoints += sequence::prefixSum(flag3, 0, nn);
    if (numCells > cellCapacity) {
      cout << "error, grid insert exceeded cell capacity, abort()" << endl;abort();}
    if(freeFlag) free(flag);
    if(freeFlag3) free(flag3);
  }

  //return no-duplicate non-sparse points in nghoods of PP
  _seq<objT> checkConflictSerial(objT* PP, intT nn, objT* down=NULL) {
    if(!down) down = newA(objT, nghSizeSerial(PP, nn));
    auto sparseFilter = [&](objT p) {
                          if(!isSparse(p)) return true;
                          else return false;};
    auto fGet = [&](intT i, objT p) {
                  down[i] = p;};
    intT total = nghPointMapSerial(PP, nn, sparseFilter, fGet);
    return _seq<objT>(down, total);
  }

  _seq<objT> checkConflictParallel(objT* PP, intT nn, objT* down=NULL, intT* flag=NULL) {
    if(!down) down = newA(objT, nghSizeParallel(PP, nn, flag));
    auto sparseFilter = [&](objT p) {
                          if(!isSparse(p)) return true;
                          else return false;};
    auto fGet = [&](intT i, objT p) {
                  down[i] = p;};
    intT total = nghPointMapParallel(PP, nn, sparseFilter, fGet, flag);
    //intT total = nghPointMapHashParallel(PP, nn, sparseFilter, fGet, false, flag);
    return _seq<objT>(down, total);
  }

  //return sparse points among PP
  _seq<objT> checkSparseSerial(objT* PP, intT nn, objT* P1=NULL) {
    if (nn<=0) return _seq<objT>();
    if(!P1) P1=newA(objT, nn);
    intT ii=0;
    for (intT i=0; i<nn; ++i) {
      if(isSparse(PP[i])) P1[ii++]=PP[i];}
    return _seq<objT>(P1,ii);}

  _seq<objT> checkSparseParallel(objT* PP, intT nn, intT* flag=NULL, objT* P1=NULL) {
    if (nn<=0) return _seq<objT>();
    bool freeFlag=false;
    if(!flag) {
      flag=newA(intT, nn+1);
      freeFlag=true;}
    parallel_for(0, nn, [&](intT i) {
	if(isSparse(PP[i])) flag[i]=1;
	else flag[i]=0;
      });
    intT ii = sequence::prefixSum(flag, 0, nn);
    flag[nn] = ii;
    if(!P1) P1=newA(objT, nn);
    parallel_for(0, nn, [&](intT i) {
	if(flag[i] != flag[i+1]) P1[flag[i]]=PP[i];
      });
    if(freeFlag) free(flag);
    return _seq<objT>(P1,ii);}

  /**
  *   Serially, return all points in the grid.
  *   @param PP auxiliary memory of size size().
  *   @return sequence containing all the points (using PP).
  */
  _seq<objT> collectPointsSerial(objT* PP=NULL) {
    if(!PP) PP=newA(objT, size());
    auto fGet = [&](intT i, objT p) {
                  PP[i] = p;};
    intT total = allPointMapSerial(fGet);
    return _seq<objT>(PP, total);
  }

  /**
  *   In parallel, return all points in the grid.
  *   @param PP auxiliary memory of size size().
  *   @param flag auxiliary memory of size size().
  *   @return ...
  */
  _seq<objT> collectPointsParallel(objT* PP=NULL, intT* flag=NULL) {
    bool freeFlag=false;
    if(!flag) {
      flag=newA(intT, size());
      freeFlag=true;}
    if(!PP) PP=newA(objT, size());
    auto fGet = [&](intT i, objT p) {
                  PP[i] = p;};
    intT total = allPointMapParallel(fGet);
    if(freeFlag)free(flag);
    return _seq<objT>(PP, total);
  }

  void printCells() {
    for (intT c=0; c<numCells; ++c) {
      cout << "cell " << c << ": ";
      cells[c].printPoints();
    }
  }

};

#endif
