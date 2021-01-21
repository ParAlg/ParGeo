#ifndef KDT_H
#define KDT_H

#include "parlay/parallel.h"
#include "parlay/sequence.h"

//make dim intrinsic to objt todo
template<int dim, class objT>
class kdNode {
  typedef int intT;
  typedef double floatT;
  typedef point<dim> pointT;
  typedef kdNode<dim, objT> nodeT;

  static const int boxInclude = 0;
  static const int boxOverlap = 1;
  static const int boxExclude = 2;

  // Data fields

  int k;
  pointT pMin, pMax;
  parlay::slice<objT**, objT**> items;
  intT n;
  nodeT* left;
  nodeT* right;
  nodeT* sib;

  // Methods

  inline void boundingBoxSerial() {
    pMin = pointT(items[0]->coordinate());
    pMax = pointT(items[0]->coordinate());
    for(intT i=0; i<n; ++i) {
      pMin.minCoords(items[i]->coordinate());
      pMax.maxCoords(items[i]->coordinate());
    }}

  inline void boundingBoxParallel() {
    intT P = parlay::num_workers()*8;
    intT blockSize = (n+P-1)/P;
    pointT localMin[P];
    pointT localMax[P];
    for (intT i=0; i<P; ++i) {
      localMin[i] = pointT(items[0]->coordinate());
      localMax[i] = pointT(items[0]->coordinate());}
    parlay::parallel_for(0, P,
  		 [&](intT p) {
  		   intT s = p*blockSize;
  		   intT e = min((intT)(p+1)*blockSize,n);
  		   for (intT j=s; j<e; ++j) {
  		     localMin[p].minCoords(items[j]->coordinate());
  		     localMax[p].maxCoords(items[j]->coordinate());}
  		 });
    pMin = pointT(items[0]->coordinate());
    pMax = pointT(items[0]->coordinate());
    for(intT p=0; p<P; ++p) {
      pMin.minCoords(localMin[p].x);
      pMax.maxCoords(localMax[p].x);}
  }

  inline intT splitItemSerial(floatT xM) {
    if (n < 2) {
      cout << "error, kdTree splitting singleton, abort" << endl;abort();}
    intT lPt = 0;
    intT rPt = n-1;
    while (lPt < rPt) {
      if (items[lPt]->coordinate(k)>=xM) {
        while (items[rPt]->coordinate(k)>=xM && lPt < rPt) {
          rPt--;
        }
        if (lPt < rPt) {
          swap(items[lPt], items[rPt]);
          rPt--; }
        else { break;}
      }
      lPt++;
    }
    if (items[lPt]->coordinate(k) < xM) lPt++;
    return lPt;
  }

  inline int boxCompare(pointT pMin1, pointT pMax1, pointT pMin2, pointT pMax2) {
    bool exclude = false;
    bool include = true;//1 include 2
    for(int i=0; i<dim; ++i) {
      if (pMax1[i]<pMin2[i] || pMin1[i]>pMax2[i]) exclude = true;
      if (pMax1[i]<pMax2[i] || pMin1[i]>pMin2[i]) include = false;
    }
    if (exclude) return boxExclude;
    else if (include) return boxInclude;
    else return boxOverlap;
  }

  inline bool itemInBox(pointT pMin1, pointT pMax1, objT* item) {
    for(int i=0; i<dim; ++i) {
      if (pMax1[i]<item->coordinate(i) || pMin1[i]>item->coordinate(i)) return false;
    }
    return true;
  }

  intT findWidest() {
    floatT xM = -1;
    for (int kk=0; kk<dim; ++kk) {
      if (pMax[kk]-pMin[kk]>xM) {
        xM = pMax[kk]-pMin[kk];
        k = kk;}}
    return k;
  }

  void constructSerial(nodeT *space, intT leafSize) {
    boundingBoxSerial();
    sib = NULL;
    if (n <= leafSize) {
      left = NULL; right = NULL;
    } else {
      if (!space[0].isEmpty() || !space[1].isEmpty()) {
        cout << "error, kdNode overwrite, abort" << endl;abort();}

      intT k = findWidest();
      floatT xM = (pMax[k]+pMin[k])/2;

      // Split items by xM (serial)
      intT median = splitItemSerial(xM);

      if (median == 0 || median == n) {median = ceil(n/2.0);}

      // Recursive construction
      space[0] = nodeT(items.cut(0, median), median, space+1, leafSize);
      space[2*median-1] = nodeT(items.cut(median, n), n-median, space+2*median, leafSize);
      left = space;
      right = space+2*median-1;
      left->sib = right;
      right->sib = left;
    }
  }

  void constructParallel(nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize) {
    boundingBoxParallel();

    sib = NULL;
    if (n <= leafSize) {
      left = NULL; right = NULL;
    } else {
      if (!space[0].isEmpty() || !space[1].isEmpty()) {
        cout << "error, kdNode overwrite, abort" << endl;abort();}

      intT k = findWidest();
      floatT xM = (pMax[k]+pMin[k])/2;

      // Split items by xM in dim k (parallel)
      parlay::parallel_for(0, n,
		   [&](intT i) {
		     if (items[i]->coordinate(k)<xM) flags[i]=1;
		     else flags[i] = 0;});
      auto mySplit = parlay::internal::split_two(items, flags);
      auto splited = mySplit.first;
      intT median = mySplit.second;
      parlay::parallel_for(0, n, [&](intT i) {items[i] = splited[i];}); // Copy back

      if (median == 0 || median == n) {median = (n/2.0);}

      // Recursive construction
      parlay::par_do([&](){space[0] = nodeT(items.cut(0, median), median, space+1, flags.cut(0, median), leafSize);},
		     [&](){space[2*median-1] = nodeT(items.cut(median, n), n-median, space+2*median, flags.cut(median, n), leafSize);});
      left = space;
      right = space+2*median-1;
      left->sib = right;
      right->sib = left;
    }
  }

public:

  inline nodeT* L() {return left;}

  inline nodeT* R() {return right;}

  inline nodeT* siblin() {return sib;}//todo

  inline intT size() {return n;}

  inline objT* operator[](intT i) {return items[i];}

  inline void setEmpty() {n=-1;}

  inline bool isEmpty() {return n<0;}

  inline bool isLeaf() {return !left;}//check

  kdNode(parlay::slice<objT**, objT**> itemss, intT nn, nodeT *space, parlay::slice<bool*, bool*> flags, intT leafSize=16): items(itemss), n(nn) {
    if (n>2000) constructParallel(space, flags, leafSize);
    else constructSerial(space, leafSize);
  }

  kdNode(parlay::slice<objT**, objT**> itemss, intT nn, nodeT *space, intT leafSize=16): items(itemss), n(nn) {
    constructSerial(space, leafSize);//todo get rid of intT n
  }

};

template<int dim, class objT>
kdNode<dim, objT>* buildKdt(parlay::sequence<objT> P, bool parallel=true, bool noCoarsen=false) {
  typedef kdNode<dim, objT> nodeT;

  size_t n = P.size();

  auto items = parlay::sequence<objT*>(n);
  parlay::parallel_for(0, n, [&](size_t i) {items[i]=&P[i];});
  parlay::slice<objT**, objT**> itemSlice = parlay::slice(items.begin(), items.end());

  auto root = (nodeT*) malloc(sizeof(nodeT)*(2*n-1));
  parlay::parallel_for(0, 2*n-1, [&](size_t i) {
				   root[i].setEmpty();
				 });

  if (parallel) {
    auto flags = parlay::sequence<bool>(n);
    auto flagSlice = parlay::slice(flags.begin(), flags.end());
    root[0] = nodeT(itemSlice, n, root+1, flagSlice, noCoarsen ? 1 : 16);
  } else {
    root[0] = nodeT(itemSlice, n, root+1, noCoarsen ? 1 : 16);
  }

  return root;
}

#endif
