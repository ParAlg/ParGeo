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
#include <set>
#include "geometry.h"
#include "shared.h"
#include "pbbs/sequence.h"
#include "pbbs/ndHash.h"
#include "pbbs/parallel.h"
#include "pbbs/parallel.h"
#include "pbbs/gettime.h"
#include "augPoint.h"
#include "serialHeap.h"
#include "parallelHeap.h"
using namespace std;

// *************************************************************
//    Helpers
// *************************************************************

template<class T>
struct getItem {
  T *A;
  getItem(T *AA): A(AA) {}
  T operator() (intT i){
    return A[i];}
};

// *************************************************************
//   Sparse partitions
// *************************************************************

template<int dim>
struct sieve {
  typedef sieve<dim> sieveT;
  typedef point<dim> geoPointT;
  typedef augPoint<dim> pointT;
  typedef grid<dim, pointT> gridT;
  typedef pointPair<dim> pointPairT;
  typedef pointPair<dim> pivotT;
  typedef _seq<pointT> containerT;
  typedef double floatT;
  typedef serialHeap<dim> heapSerialT;
  typedef parallelHeap<dim> heapParallelT;
  static constexpr double dMax = numeric_limits<double>::max();
  static constexpr intT levelMax = 20;
  static constexpr intT heapThresh = 1000;
  static constexpr bool keepSparseSets = false;
  static constexpr bool noRandom = false;

  struct statistics {
    double getPivot = 0;
    double gridInsert = 0;
    double gridInsertCheck = 0;
    double heapInsert = 0;
    double insertRebuild = 0;
    intT insertRebuildCount = 0;
    double totalInsert = 0;
    intT insertCalls = 0;

    double gridErase = 0;
    double gridEraseCheck = 0;
    double gridEraseMove = 0;
    double heapErase = 0;
    double eraseRebuild = 0;
    intT eraseRebuildCount = 0;
    double totalErase = 0;
    intT eraseCalls = 0;
  };
  struct statistics stats;
  void printStats() {
    std::cout << std::setprecision(3);
    cout << ">>> dynamic calls stats" << endl << endl;
    cout << "total insert calls = " << stats.insertCalls << endl;
    cout << "total insert time = " << stats.totalInsert << endl;
    cout << "  find pivot time = " << stats.getPivot;
    cout << " (" << intT(100*stats.getPivot/stats.totalInsert) << "%)" << endl;
    cout << "  grid insert time = " << stats.gridInsert;
    cout << " (" << intT(100*stats.gridInsert/stats.totalInsert) << "%)" << endl;
    cout << "  grid check time = " << stats.gridInsertCheck;
    cout << " (" << intT(100*stats.gridInsertCheck/stats.totalInsert) << "%)" << endl;
    cout << "  heap time = " << stats.heapInsert;
    cout << " (" << intT(100*stats.heapInsert/stats.totalInsert) << "%)" << endl;
    cout << "  rebuild time = " << stats.insertRebuild;
    cout << " (" << intT(100*stats.insertRebuild/stats.totalInsert) << "%)";
    cout << " (counts = " << stats.insertRebuildCount << ")" << endl;
    cout << endl;
    cout << "total erase calls = " << stats.eraseCalls << endl;
    cout << "total erase time = " << stats.totalErase << endl;
    cout << "  grid erase time = " << stats.gridErase;
    cout << " (" << intT(100*stats.gridErase/stats.totalErase) << "%)" << endl;
    cout << "  grid check time = " << stats.gridEraseCheck;
    cout << " (" << intT(100*stats.gridEraseCheck/stats.totalErase) << "%)" << endl;
    cout << "  grid move time = " << stats.gridEraseMove;
    cout << " (" << intT(100*stats.gridEraseMove/stats.totalErase) << "%)" << endl;
    cout << "  heap time = " << stats.heapErase;
    cout << " (" << intT(100*stats.heapErase/stats.totalErase) << "%)" << endl;
    cout << "  rebuild time = " << stats.eraseRebuild;
    cout << " (" << intT(100*stats.eraseRebuild/stats.totalErase) << "%)";
    cout << " (counts = " << stats.eraseRebuildCount << ")" << endl;
    cout << endl;
    std::cout << std::setprecision(6);
  }

  intT n;
  gridT** fullSets=NULL;
  gridT** sparseSets=NULL;
  pivotT* pivots=NULL;
  intT numLevels;
  intT cellMax;
  geoPointT pMin;
  sieve(geoPointT pMinn, intT maxN): pMin(pMinn), cellMax(maxN), numLevels(0) {
    fullSets = newA(gridT*, levelMax*2+1);
    for(intT i=0; i<levelMax*2+1; ++i) fullSets[i]=NULL;
    sparseSets = fullSets+levelMax;
    refreshGP(maxN);//todo, more efficient
    refreshGFlag(maxN);//todo, more efficient
    pivots = newA(pivotT, levelMax+1);
    for(intT i=0; i<levelMax+1; ++i) pivots[i]=pivotT();
  }
  ~sieve() {
    if(fullSets) free(fullSets);//sparseSets included
  }

  intT rands[10] = {1,111,1100,1012,300123,1230123,12031230,1231239,5430,12093812};

  inline pivotT getPivot(intT l) {//todo in progress
    return pivots[l];
  }

  inline void assignPivot(intT l, pivotT vv) {
    pivots[l] = vv;
  }

  //get pivot without using grid (random selection + linear scan)
  pivotT getPivotSerial(pointT* PP, intT nn) {
    static const intT numTry = 10;
    pivotT results[numTry];
    intT choices[numTry];
    for(intT i=0; i<numTry; ++i) {
      if(noRandom) choices[i] = rands[i]%nn;
      else choices[i] = rand()%nn;
      results[i] = pivotT();}
    for(intT j=0; j<numTry; ++j) {
      for(intT i=0; i<nn; ++i) {
        auto tmp = pivotT(PP[choices[j]].pt(), PP[i].pt());
        if (tmp.dist >0) results[j].closer(tmp);}
    }
    auto result = pivotT();
    result.dist = 0;
    for(intT i=0; i<numTry; ++i) {
      if (!results[i].isEmpty()) {
        result = result.dist > results[i].dist ? result : results[i];}
    }
    if (result.isEmpty() || result.dist == 0) {
      cout << "nn = " << nn << endl;
      cout << "error, getPivotSerial, empty pivot, abort" << endl; abort();
    } else return result;}

  pivotT getPivotParallel(pointT* PP, intT nn) {
    static const intT numTry = 10;
    pivotT results[numTry];
    intT choices[numTry];
    for(intT i=0; i<numTry; ++i) {
      if(noRandom) choices[i] = rands[i]%nn;
      else choices[i] = rand()%nn;
      results[i] = pivotT();}
    floatT* xDists = newA(floatT, nn);
    for (intT j=0; j<numTry; ++j) {
      par_for(intT i=0; i<nn; ++i) {
        xDists[i] = PP[i].pointDist(PP[choices[j]]);
      }
      xDists[choices[j]] = dMax;
      intT iMin = sequence::minIndex<floatT,intT,getItem<floatT>>(0, nn, getItem<floatT>(xDists));
      results[j] = pivotT(PP[choices[j]].pt(), PP[iMin].pt(), xDists[iMin]);
    }
    auto result = pivotT();
    result.dist = 0;
    for(intT i=0; i<numTry; ++i) {
      if (!results[i].isEmpty()) {
        result = result.dist > results[i].dist ? result : results[i];}
    }
    free(xDists);
    if (result.isEmpty() || result.dist==0) {
      cout << "nn = " << nn << endl;
      cout << "error, getPivotParallel, zero, abort" << endl; abort();
    } else return result;}

  inline void newFullSet(intT l, pivotT myPivot) {
    assignPivot(l, myPivot);
    fullSets[l] = new gridT(cellMax, pMin, myPivot.dist/(6*dim));
  }

  inline void newSparseSet(intT l, pivotT myPivot) {
    assignPivot(l, myPivot);//should be consistent with fullset
    sparseSets[l] = new gridT(cellMax, pMin, myPivot.dist/(6*dim));
  }

  intT *gFlag = NULL;
  intT gFlagSize = -1;
  void refreshGFlag(intT nn) {
    if (!gFlag || gFlagSize < nn) {
      if(gFlag) free(gFlag);
      gFlagSize = nn*4+4;
      gFlag=newA(intT, gFlagSize);}
  }
  inline intT* gFlag1() {return gFlag;};
  inline intT* gFlag2() {return gFlag+gFlagSize/2;};

  pointT *gP = NULL;
  intT gPSize = -1;
  void refreshGP(intT nn) {
    if (!gP || gPSize < nn) {
      if(gP) free(gP);
      gPSize = nn*2;
      gP=newA(pointT, gPSize);}
  }
  inline pointT* gP1() {return gP;}
  inline pointT* gP2() {return gP+gPSize/2;}

  //find closest pair of a point in partition
  inline pointPairT findPair(gridT* g, pointT p) {
    auto R = pointPairT();
    auto fPair = [&] (pointT& q) {
                   if(p!=q) R.closer(p.pt(),q.pt());
                   return false;
                 };
    g->nghPointMap(p.coordinate(), fPair);
    if (R.isEmpty()) R.u=p.pt();
    return R;
  }

  //get no-duplicate neighbors for input point array
  _seq<geoPointT> getNeighborsParallel(gridT* g, pointT* PP, intT nn, pointT* P2=NULL, intT* flag=NULL) {
    if(nn==0) return _seq<geoPointT>();
    bool freeP2=false;
    if(!P2) {
      P2 = newA(pointT, g->nghSizeParallel(PP, nn, flag));
      freeP2=true;}
    auto fGet = [&](intT i, pointT p) {
                     P2[i] = p;};
    intT total = g->nghPointMapParallel(PP, nn, [&](pointT p){return true;}, fGet, flag);
    //intT total = g->nghPointMapHashParallel(PP, nn, [&](pointT p){return true;}, fGet, false, flag);
    auto A = newA(geoPointT, total);
    par_for(intT i=0; i<total; ++i) {A[i] = P2[i].pt();}
    if(freeP2) free(P2);
    return _seq<geoPointT>(A, total);
  }

  //call after full set insert
  bool needInsertionRebuild(intT l, gridT* g, bool serial=false) {
    intT numPoints=g->size();

    if (!noRandom && rand() % (numPoints) == 0) return true;

    //auto oldPivot = g->getPivot();
    auto oldPivot = getPivot(l);
    auto newPair = findPair(g, oldPivot.u);
    if (newPair.dist < oldPivot.dist) return true;
    else return false;
  }

  bool insert(pointT* PP, intT nn, bool serial, intT l=0, bool rebuild=true, bool recordStats=true) {
    static const bool verbose = false;

    if (nn==0) return false;

    timing t0, tt;
    if(recordStats) {
      stats.insertCalls++;
      tt.start();}

    pointT* pInsert = PP;
    intT nInsert = nn;
    bool everRebuilt = false;
    bool grew = false;
    do {
      if (verbose) cout << "-- insert to level " << l << ", " << nInsert << " points" << endl;
      if (fullSets[l] == NULL) {
        if (verbose) cout << "create new grid" << endl;
        pivotT pivot;

        if(recordStats) t0.start();
        if(serial) pivot = getPivotSerial(pInsert, nInsert);
        else pivot = getPivotParallel(pInsert, nInsert);
        if(recordStats) stats.getPivot += t0.stop();

        newFullSet(l, pivot);
        if (keepSparseSets) newSparseSet(l, pivot);
        numLevels++;
        grew = true;
      }
      if (verbose) printLevel(l);
      bool roundRebuilt=false;
      bool needRebuild=false;

      if(recordStats) t0.start();
      if (serial) {
        fullSets[l]->insertSerial(pInsert, nInsert);
        if(rebuild) needRebuild = needInsertionRebuild(l, fullSets[l], true);
      } else {
        fullSets[l]->insertParallel(pInsert, nInsert, gFlag1(), gFlag2());
        if(rebuild) needRebuild = needInsertionRebuild(l, fullSets[l], false);}
      if(recordStats) stats.gridInsert += t0.stop();

      if (rebuild && needRebuild) {
        if(recordStats) t0.start();
        if(recordStats) stats.insertRebuildCount += 1;
        roundRebuilt=true;
        everRebuilt = true;
        if (verbose) cout << "rebuild at lvl." << l << endl;
        containerT collected;
        if(serial) collected = fullSets[l]->collectPointsSerial(gP2());//must use gP2
        else collected = fullSets[l]->collectPointsParallel(gP2(), gFlag1());//must use gP2
        pInsert = collected.A; nInsert = collected.n;
        intT ll = 0;
        while(fullSets[l+ll] != NULL) {//delete levels >=l
          if(verbose)cout << "deleting grid " << l+ll << endl;
          delete fullSets[l+ll];
          fullSets[l+ll] = NULL;
          if(keepSparseSets) {
            delete sparseSets[l+ll];
            sparseSets[l+ll] = NULL;}
          ll++;
          numLevels--;}
        deleteHeap();
        insert(pInsert, nInsert, serial, l, rebuild, true);//rebuild from l
        if(recordStats) stats.insertRebuild += t0.stop();
        return true;
      }

      if(keepSparseSets) {
        if(recordStats) t0.start();
        auto candidates = fullSets[l]->collectPointsSerial();
        auto sps = fullSets[l]->checkSparseSerial(candidates.A, candidates.n);
        if (serial) sparseSets[l]->insertSerial(sps.A, sps.n);
        else sparseSets[l]->insertParallel(sps.A, sps.n, gFlag1(), gFlag2());
        free(candidates.A); free(sps.A);
        if(recordStats) stats.gridInsert += t0.stop();
      }

      if(recordStats) t0.start();
      if (l == heapLevel) {
        if(verbose) cout << "heap insert, size before = " << heapSize() << endl;
        if(serial) heapInsertSerial(pInsert, nInsert);
        else heapInsertParallel(pInsert, nInsert, gP1(), gFlag1());
        if(verbose) cout << " size after = " << heapSize() << endl;
      }
      if(recordStats) stats.heapInsert += t0.stop();

      if(recordStats) t0.start();
      containerT down;
      if (serial) down = fullSets[l]->checkConflictSerial(pInsert, nInsert, gP1());
      else down = fullSets[l]->checkConflictParallel(pInsert, nInsert, gP1(), gFlag1());
      if(recordStats) stats.gridInsertCheck += t0.stop();

      if (verbose) cout << "|down_i| = " << down.n << endl;
      if (verbose) printLevel(l);
      if (verbose) cout << endl;

      pInsert = down.A; nInsert = down.n;
      l++;
      numLevels = max(numLevels, l);
    } while (l<levelMax && nInsert > 0);

    if(recordStats) t0.start();
    if((!hasHeap() && numLevels>0) || // no heap yet
       (grew && (hasHeap()&&heapSize()>heapThresh) && whichPartition(lastPartition())>heapLevel)) {
      if(verbose) cout << "rebuild heap first time or for optimization" << endl << endl;
      if(serial) reheapSerial(gP1(), gP2());
      else reheapParallel(gP1(), gP2(), gFlag1(), gFlag2());
    }
    if(recordStats) stats.heapInsert += t0.stop();

    if(recordStats) stats.totalInsert += tt.stop();
    return everRebuilt;
  }

  bool needEraseRebuildSerial(pivotT pivot, pointT* PP, intT nn) {
    bool need = false;
    for(intT i=0; i<nn; ++i) {
      if (samePoint<dim>(PP[i].x(), pivot.u.x)) {
        need=true;
        break;}
    }
    return need;}

  bool needEraseRebuildParallel(pivotT pivot, pointT* PP, intT nn) {
    bool need = false;
    par_for(intT i=0; i<nn; ++i) {
      if (samePoint<dim>(PP[i].x(), pivot.u.x)) need=true;
    }
    return need;}

  bool erase(pointT* PP, intT nn, bool serial=false) {
    static const bool rebuild = true;
    static const bool verbose = false;

    if (nn==0) return false;

    stats.eraseCalls += 1;
    timing t0, tt;
    tt.start();

    intT l=numLevels-1;
    bool everRebuilt = false;
    containerT up;
    bool needRebuild1, needRebuild2=false;
    intT rebuildAtLevel = -1;
    do {
      if (l>=0) {
        if (verbose) cout << "-- erase from level " << l << endl;
        if (verbose) cout << "nn = " << nn << endl;
        if (verbose) printLevel(l);

        t0.start();
        containerT upNeighbors;
        if (serial) {
          upNeighbors = fullSets[l]->eraseSerial(PP, nn, true, gP1());
          if(keepSparseSets) sparseSets[l]->eraseSerial(PP, nn, false);
        } else {
          upNeighbors = fullSets[l]->eraseParallel(PP, nn, gFlag1(), true, gP1());
          if(keepSparseSets) sparseSets[l]->eraseParallel(PP, nn, gFlag1(), false);}
        stats.gridErase += t0.stop();

        t0.start();
        if(serial) up = fullSets[l]->checkSparseSerial(upNeighbors.A, upNeighbors.n, gP2());
        else up = fullSets[l]->checkSparseParallel(upNeighbors.A, upNeighbors.n, gFlag1(), gP2());
        stats.gridEraseCheck += t0.stop();

        t0.start();
        if(l==heapLevel) {
          if(serial) heapEraseSerial(PP, nn);
          else heapEraseParallel(PP, nn, gP1(), gFlag1());
        }
        stats.heapErase += t0.stop();

        t0.start();
        needRebuild1 = false;
        if (rebuild && serial && needEraseRebuildSerial(getPivot(l), PP, nn)) {
          needRebuild1=true;
        } else if (rebuild && !serial && needEraseRebuildParallel(getPivot(l), PP, nn)) {
          needRebuild1=true;}
        stats.gridErase += t0.stop();

        if (l<numLevels-1) {
          if(verbose) cout << "move |up|=" << up.n << " from " << l+1 << " to " << l  << endl;
          if(verbose) printLevel(l+1);

          t0.start();
          if (serial) {
            fullSets[l+1]->eraseSerial(up.A, up.n, false);
            if(keepSparseSets) {
              sparseSets[l+1]->eraseSerial(up.A, up.n, false);
              sparseSets[l]->insertSerial(up.A, up.n);}
          } else {
            fullSets[l+1]->eraseParallel(up.A, up.n, gFlag1(), false);
            if(keepSparseSets) {
              sparseSets[l+1]->eraseParallel(up.A, up.n, gFlag1(), false);
              sparseSets[l]->insertParallel(up.A, up.n, gFlag1(), gFlag2());}
          }
          stats.gridEraseMove += t0.stop();

          t0.start();
          if (l==heapLevel) {
            if(serial) heapInsertSerial(up.A, up.n);
            else heapInsertParallel(up.A, up.n, gP1(), gFlag1());
          }
          if (l+1==heapLevel) {
            if(serial) heapEraseSerial(up.A, up.n);
            else heapEraseParallel(up.A, up.n, gP1(), gFlag1());
          }
          stats.heapErase += t0.stop();
        }
      }

      if(l<numLevels-1) {
        t0.start();
        if (rebuild && serial && needEraseRebuildSerial(getPivot(l+1), up.A, up.n)) {
          needRebuild2 |= true;
        } else if (rebuild && !serial && needEraseRebuildParallel(getPivot(l+1), up.A, up.n)) {
          needRebuild2 |= true;}
        stats.gridErase += t0.stop();

        if (rebuild && needRebuild2) rebuildAtLevel = l+1;

        if (verbose) printLevel(l+1);
      }

      if (verbose) cout << endl;
      needRebuild2 = needRebuild1;
      l--;
    } while(l>=-1);

    bool shrinked = false;
    if (rebuild && rebuildAtLevel >= 0) {
      t0.start();
      stats.eraseRebuildCount += 1;
      if(verbose) cout << "rebuild from lvl. " << rebuildAtLevel << " due to erase" << endl;
      containerT collected;
      if(serial) collected = fullSets[rebuildAtLevel]->collectPointsSerial(gP2());
      else collected = fullSets[rebuildAtLevel]->collectPointsParallel(gP2(), gFlag1());
      intT ll = 0;
      while(fullSets[rebuildAtLevel+ll] != NULL) {
        delete fullSets[rebuildAtLevel+ll];
        fullSets[rebuildAtLevel+ll] = NULL;
        if(keepSparseSets) {
          delete sparseSets[rebuildAtLevel+ll];
          sparseSets[rebuildAtLevel+ll] = NULL;}
        ll++;
        numLevels--;
      }
      if(verbose) cout << "rebuild size = " << collected.n << " at level " << rebuildAtLevel << endl;
      deleteHeap();
      insert(collected.A, collected.n, serial, rebuildAtLevel, false, false);
      stats.eraseRebuild += t0.stop();
    }

    if(rebuildAtLevel < 0) {
      //if never rebuilt, check if had fewer levels
      intT ll=numLevels-1;
      while (ll>=0) {
        intT nextSize = fullSets[ll]->size();
        if (nextSize<=0) {//level removed
          delete fullSets[ll]; fullSets[ll] = NULL;
          if (keepSparseSets) {
            delete sparseSets[ll]; sparseSets[ll] = NULL;}
          numLevels--;
          shrinked=true;}
        ll--;
      }
    }

    t0.start();
    if(rebuildAtLevel>=0 ||
       (shrinked && whichPartition(lastPartition())<heapLevel && numLevels>0)) {
      if(verbose) cout << "rebuild heap due to erase" << endl << endl;
      if(serial) reheapSerial(gP1(), gP2());
      else reheapParallel(gP1(), gP2(), gFlag1(), gFlag2());
    }
    stats.heapErase += t0.stop();

    stats.totalErase += tt.stop();
    return false;
  }

  heapSerialT *heapSerial=NULL;
  heapParallelT *heapParallel=NULL;
  intT heapLevel = -1;

  void heapStats() {
    if(heapSerial) return;
    else if(heapParallel) heapParallel->printStats();;
  }

  inline intT heapSize() {
    if(heapSerial) return heapSerial->size();
    else if(heapParallel) return heapParallel->size();
    else return -1;
  }

  inline pointPairT heapMin() {
    if(heapSerial) return heapSerial->getMin();
    if(heapParallel) return heapParallel->getMin();
    return pointPairT();
  }

  inline void deleteHeap() {
    if(heapSerial) {
      delete heapSerial;
      heapSerial = NULL;}
    if(heapParallel) {
      delete heapParallel;
      heapParallel = NULL;}
  }

  inline bool hasHeap() {
    return heapSerial || heapParallel;
  }

  void heapInsertSerial(pointT* PP, intT nn) {
    auto g = fullSets[heapLevel];
    if (g==NULL || heapSerial==NULL || nn<=0) return;
    for (intT i=0; i<nn; ++i) {
      pointT u = PP[i];
      auto cp = findPair(g, u.pt());
      heapSerial->insert(cp);
      intT repeats = 0;
      auto fUpdate = [&] (pointT& v) {
                    if (repeats<=0 && v==u) {
                      return false;repeats++;}
                    heapSerial->decreaseKey(pointPairT(v.pt(),u.pt()));
                    return false;};
      g->nghPointMap(u.coordinate(), fUpdate);
    }
  }

  //todo, create serial routine for smaller inputs
  void heapInsertParallel(pointT* PP, intT nn, pointT* P2=NULL, intT* flag=NULL) {
    auto g = fullSets[heapLevel];
    if (g==NULL || heapParallel==NULL || nn<=0) return;
    _seq<geoPointT> neighbors = getNeighborsParallel(g, PP, nn, P2, flag);
    heapParallel->erase(neighbors.A, neighbors.n);
    auto RR = newA(pointPairT, neighbors.n);
    par_for(intT i=0; i<neighbors.n; ++i) {
      RR[i] = findPair(g, neighbors.A[i]);
    }
    heapParallel->insert(RR, neighbors.n);
  }

  void heapEraseSerial(pointT* PP, intT nn) {
    auto g = fullSets[heapLevel];
    if (g==NULL || heapSerial==NULL) return;
    for (intT i=0; i<nn; ++i) {
      pointT u = PP[i];
      heapSerial->erase(u.pt());
      auto fReplace = [&] (pointT& v) {
                        floatT oldDist = heapSerial->dist(v.pt());
                        auto Rn = findPair(g, v.pt());
                        if (Rn.dist > oldDist) heapSerial->replace(Rn);
                        return false;};
      g->nghPointMap(u.coordinate(), fReplace);
    }
  }

  void heapEraseParallel(pointT* PP, intT nn, pointT* P2=NULL, intT* flag=NULL) {
    auto g = fullSets[heapLevel];
    if (g==NULL || heapParallel==NULL) return;
    auto A = newA(geoPointT, nn);//todo
    par_for(intT i=0; i<nn; ++i) {
      A[i]=PP[i].pt();
    }
    heapParallel->erase(A, nn);
    free(A);
    _seq<geoPointT> neighbors = getNeighborsParallel(g, PP, nn, P2, flag);
    heapParallel->erase(neighbors.A, neighbors.n);
    auto RR = newA(pointPairT, neighbors.n);
    par_for(intT i=0; i<neighbors.n; ++i) {
      RR[i] = findPair(g, neighbors.A[i]);
    }
    heapParallel->insert(RR, neighbors.n);
  }

  void reheapSerial(pointT* gP1, pointT* gP2) {
    deleteHeap();
    heapLevel = whichPartition(lastPartition());
    if (!fullSets[heapLevel] || numLevels<=1) return;
    intT nn = fullSets[heapLevel]->size();
    auto Rs = newA(pointPairT, nn);
    auto fCp = [&](intT i, pointT p) {
                 Rs[i] = findPair(fullSets[heapLevel], p);};
    fullSets[heapLevel]->allPointMapSerial(fCp);

    heapSerial = new heapSerialT(cellMax, fullSets[heapLevel]->pMin, fullSets[heapLevel]->r);
    for(intT i=0; i<nn; ++i) {
      heapSerial->insert(Rs[i]);}
    free(Rs);
  }

  void reheapParallel(pointT* gP1, pointT* gP2, intT* gFlag1, intT* gFlag2) {
    heapLevel = whichPartition(lastPartition());
    if (!fullSets[heapLevel] || numLevels<=1) return;
    intT nn = fullSets[heapLevel]->size();
    auto Rs = newA(pointPairT, nn);
    auto fCp = [&](intT i, pointT p) {
                 Rs[i] = findPair(fullSets[heapLevel], p);};
    fullSets[heapLevel]->allPointMapParallel(fCp, gFlag1);
    if(heapParallel) {
      auto stats = heapParallel->stats;
      deleteHeap();
      heapParallel = new heapParallelT(cellMax, fullSets[heapLevel]->pMin, fullSets[heapLevel]->r);
      heapParallel->stats = stats;
    } else {
      heapParallel = new heapParallelT(cellMax, fullSets[heapLevel]->pMin, fullSets[heapLevel]->r);
    }
    heapParallel->insert(Rs, nn);
    free(Rs);
  }

  pointPairT closestPairGridSerial(intT use) {
    return fullSets[heapLevel]->closestPairSerial();}

  pointPairT closestPairGridParallel(intT use) {
    return fullSets[heapLevel]->closestPairParallel();}

  pointPairT closestPairHeap() {
    if (!hasHeap()) {
      cout << "query heap = NULL" << endl;
      return pointPairT();
    } else {
      // cout << "query heap size = " << heapSize() << ", heap level = " << heapLevel << endl;
      // if (heapSize() > 0 && heapSize() != fullSets[heapLevel]->size()) {//todo, remove
      //   cout << "heap size wrong" << endl;
      //   abort();
      // }
      return heapMin();}
  }

  inline intT lastPartition() {
    intT last=numLevels-1;
    for (;last>=0;last--) {
      if (fullSets[last]->size() > 0) {
        break;}
    }
    return last;}

  inline intT whichPartition(intT last) {
    intT use = max((intT)0,last-(intT)ceil(log(2*sqrt(dim))/log(3)));
    return use;}

  pointPairT closestPairSerial() {
    return closestPairHeap();
  }

  pointPairT closestPairParallel() {
    return closestPairHeap();
  }

  void printLevel(intT l) {
    if (fullSets[l]) {
      cout << "r = " << fullSets[l]->r << endl;
      cout << "|S_" << l << "| = " << fullSets[l]->size() << endl;
      if (keepSparseSets) cout << "|S_" << l << "'| = " << sparseSets[l]->size() << endl;
    } else {
      cout << "|S_" << l << "| = NULL" << endl;
      if (keepSparseSets) cout << "|S_" << l << "'| = NULL" << endl;
    }
  }

  void printSummary() {
    std::cout << std::setprecision(2);
    cout << "nested partitions = " << endl;
    for (intT i=0; i<numLevels; ++i) {
      cout << " lvl." << i << " ";
      cout << "r=" << fullSets[i]->r << " ";
      cout << "#cells=" << fullSets[i]->numCells << " ";
      cout << "#pts=" << fullSets[i]->size();
      cout << endl;
    }
    std::cout << std::setprecision(6);
  }

};

// *************************************************************
//    DRIVER
// *************************************************************
template<int dim>
pair<point<dim>, point<dim>> closestPair(point<dim>* P, intT n, intT batches) {
  static const bool verbose = true;
  static const bool serial = false;
  cout << "dynamic closestPair of " << n << ", dim " << dim << " points" << endl << endl;
  if (n < 2) abort();

  typedef double floatT;
  typedef pointPair<dim> pointPairT;
  typedef augPoint<dim> pointT;
  typedef point<dim> geoPointT;

  geoPointT pMin;
  if (serial) pMin = pMinSerial(P, n);
  else pMin = pMinParallel(P, n);

  auto PP = newA(pointT, n);
  par_for(intT i=0; i<n; ++i) {
    PP[i] = pointT(P[i]);
  }

  intT batchSize = n/batches;
  cout << "dynamic insert+delete " << batches << " batches of " << batchSize << endl << endl;

  timing t0; t0.start();
  double tInsert = 0;
  double tErase = 0;
  double tQuery = 0;
  double tmp;
  auto sp = sieve<dim>(pMin, n);//adjust cellMax

  intT inserted = 0;
  intT queries = 0;
  for (intT b=0; b<batches; ++b) {
    if(verbose) cout << "-- batch " << b << endl;
    if (b==batches-1) {
      sp.insert(PP+inserted, max(batchSize, n-inserted), serial);
      inserted += max(batchSize, n-inserted);
    } else {
      static const bool consistencyChecker=false;

      if (!consistencyChecker) {
        sp.insert(PP+inserted, batchSize, serial);
      } else {
        intT checks = 4;
        intT vecLen = 6;
        intT fullSizes[checks*vecLen];
        intT fullCounter[checks];
        intT sparseSizes[checks*vecLen];
        floatT floatSizes[checks*vecLen];
        bool hasRebuild = false;
        for (intT i=0; i<checks; ++i) {
          fullCounter[i]=0;
          for (intT j=0; j<vecLen; ++j) {
            floatSizes[i*vecLen+j]=0;
            fullSizes[i*vecLen+j]=0;
            sparseSizes[i*vecLen+j]=0;}
        }

        intT check = 0;
        for(intT i=0; i<sp.numLevels; ++i) {
          floatSizes[check*vecLen+i]=sp.fullSets[i]->r;
          fullSizes[check*vecLen+i] = sp.fullSets[i]->size();
          if(sp.keepSparseSets)sparseSizes[check*vecLen+i] = sp.sparseSets[i]->size();}
        fullCounter[check] = sp.numLevels;

        hasRebuild |= sp.insert(PP+inserted, batchSize, serial);
        check = 1;
        for(intT i=0; i<sp.numLevels; ++i) {
          floatSizes[check*vecLen+i]=sp.fullSets[i]->r;
          fullSizes[check*vecLen+i] = sp.fullSets[i]->size();
          if(sp.keepSparseSets)sparseSizes[check*vecLen+i] = sp.sparseSets[i]->size();}
        fullCounter[check] = sp.numLevels;

        hasRebuild |= sp.erase(PP+inserted, batchSize, serial);
        check = 2;
        for(intT i=0; i<sp.numLevels; ++i) {
          floatSizes[check*vecLen+i]=sp.fullSets[i]->r;
          fullSizes[check*vecLen+i] = sp.fullSets[i]->size();
          if(sp.keepSparseSets)sparseSizes[check*vecLen+i] = sp.sparseSets[i]->size();}
        fullCounter[check] = sp.numLevels;

        hasRebuild |= sp.insert(PP+inserted, batchSize, serial);
        check = 3;
        for(intT i=0; i<sp.numLevels; ++i) {
          floatSizes[check*vecLen+i]=sp.fullSets[i]->r;
          fullSizes[check*vecLen+i] = sp.fullSets[i]->size();
          if(sp.keepSparseSets)sparseSizes[check*vecLen+i] = sp.sparseSets[i]->size();}
        fullCounter[check] = sp.numLevels;

        //printout
        cout << ">>>>>>>>>>>>>>>" << endl;
        check=0;
        for (intT i=0; i<fullCounter[check]; ++i) {
          cout << "level " << i << ", r = " << floatSizes[check*vecLen+i] <<  ": ";
          if(sp.keepSparseSets)cout << sparseSizes[check*vecLen+i] << "/";
          cout << fullSizes[check*vecLen+i] << endl;}
        cout << ">>insert " << batchSize << endl;
        check=1;
        for (intT i=0; i<fullCounter[check]; ++i) {
          cout << "level " << i << ", r = " << floatSizes[check*vecLen+i] <<  ": ";
          if(sp.keepSparseSets)cout << sparseSizes[check*vecLen+i] << "/";
          cout << fullSizes[check*vecLen+i] << endl;}
        cout << ">>erase " << batchSize << endl;
        check=2;
        for (intT i=0; i<fullCounter[check]; ++i) {
          cout << "level " << i << ", r = " << floatSizes[check*vecLen+i] <<  ": ";
          if(sp.keepSparseSets)cout << sparseSizes[check*vecLen+i] << "/";
          cout << fullSizes[check*vecLen+i] << endl;}
        cout << ">>insert " << batchSize << endl;
        check=3;
        for (intT i=0; i<fullCounter[check]; ++i) {
          cout << "level " << i << ", r = " << floatSizes[check*vecLen+i] <<  ": ";
          if(sp.keepSparseSets)cout << sparseSizes[check*vecLen+i] << "/";
          cout << fullSizes[check*vecLen+i] << endl;}
        cout << ">>>>>>>>>>>>>>>" << endl;

        //check correctness, compare 0 and 2
        bool hasError=false;
        if (fullCounter[0] != fullCounter[2]) {
          cout << "insert->erase level diff, " << fullCounter[0] << " vs " << fullCounter[2] << endl;
        } else {
          for (intT i=0; i<min(fullCounter[0],fullCounter[2]); ++i) {
            if (sp.keepSparseSets&&sparseSizes[0*vecLen+i] != sparseSizes[2*vecLen+i]) hasError=true;
            if (fullSizes[0*vecLen+i] != fullSizes[2*vecLen+i]) hasError=true;
          }
        }
        if (fullCounter[1] != fullCounter[3]) {
          cout << "erase->insert level diff, " << fullCounter[1] << " vs " << fullCounter[3] << endl;
        } else {
          for (intT i=0; i<min(fullCounter[1],fullCounter[3]); ++i) {
            if (sp.keepSparseSets&&sparseSizes[1*vecLen+i] != sparseSizes[3*vecLen+i]) hasError=true;
            if (fullSizes[1*vecLen+i] != fullSizes[3*vecLen+i]) hasError=true;
          }
        }
        if (hasError) {
          cout << "inconsistency detected" << endl;
          if(hasRebuild) {
            cout << "has rebuild, continue" << endl;
          } else if(floatSizes[1*vecLen+0] != floatSizes[3*vecLen+0]) {
            cout << "pivot has changed, continue" << endl;
          } else {
            cout << "if grid->r has changed, this is normal, just run again" << endl;
            abort();
          }
        }
      }
      inserted += batchSize;}
    if(verbose) sp.printSummary();
    tmp = tInsert; tInsert += t0.next();
    pointPairT R;
    if(serial) R = sp.closestPairSerial();
    else R = sp.closestPairParallel();
    queries ++;
    tmp = tQuery; tQuery += t0.next();
    if(verbose) cout << "inserted = " << inserted << "/" << n << endl;
    if(verbose) cout << "closest-dist = " << R.dist << endl << endl;
  }

  intT erased = 0;
  for (intT b=0; b<batches; ++b) {
    if(verbose) cout << "-- batch " << b+batches << endl;//cummulative
    if (b==batches-1) {
      auto lastBatchSize = max(batchSize, n-erased);
      sp.erase(PP, lastBatchSize, serial);
      erased += lastBatchSize;
    } else {
      sp.erase(PP+n-erased-batchSize, batchSize, serial);
      erased += batchSize;
    }
    if(verbose) sp.printSummary();
    tmp = tErase; tErase += t0.next();
    pointPairT R;
    if(serial) R = sp.closestPairSerial();
    else R = sp.closestPairParallel();
    queries ++;
    tmp = tQuery; tQuery += t0.next();
    if(verbose) cout << "erased = " << erased << "/" << n << endl;
    if(verbose) cout << "closest-dist = " << R.dist << endl << endl;
  }

  cout << "total-dynamic-insert-time = " << tInsert << endl;
  cout << "total-dynamic-erase-time = " << tErase << endl;
  cout << "total-query-time = " << tQuery << endl << endl;
  cout << "avg-dynamic-insert-time = " << tInsert/inserted << endl;
  cout << "avg-dynamic-erase-time = " << tErase/erased << endl;
  cout << "avg-query-time = " << tQuery/queries << endl;
  cout << "queries = " << queries << endl << endl;

  sp.printStats();
  sp.heapStats();

  pointPairT R;
  if(serial) R = sp.closestPairSerial();
  else R = sp.closestPairParallel();
  cout << "dynamicSieve " << R.u << ", " << R.v << ", dist " << R.dist << endl << endl;

  free(PP);
  return make_pair(R.u, R.v);
}

template pair<point<2>, point<2>> closestPair<2>(point<2>*, intT, intT);
template pair<point<3>, point<3>> closestPair<3>(point<3>*, intT, intT);
template pair<point<4>, point<4>> closestPair<4>(point<4>*, intT, intT);
template pair<point<5>, point<5>> closestPair<5>(point<5>*, intT, intT);
template pair<point<6>, point<6>> closestPair<6>(point<6>*, intT, intT);
template pair<point<7>, point<7>> closestPair<7>(point<7>*, intT, intT);
template pair<point<8>, point<8>> closestPair<8>(point<8>*, intT, intT);
template pair<point<9>, point<9>> closestPair<9>(point<9>*, intT, intT);
