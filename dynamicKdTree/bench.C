#include "dynamicKdTree.h"
#include "kdTree.h"
#include "pbbs/gettime.h"
using namespace std;

// *************************************************************
//    Tree Checker
// *************************************************************

template<int dim>
inline pair<point<dim>, point<dim>> boundingBox(point<dim> **items, intT n) {
  typedef point<dim> pointT;
  auto pMin = pointT(items[0]->coordinate());
  auto pMax = pointT(items[0]->coordinate());
  for(intT i=0; i<n; ++i) {
    pMin.minCoords(items[i]->coordinate());
    pMax.maxCoords(items[i]->coordinate());
  }
  return make_pair(pMin, pMax);
}

template<class nodeT, int dim>
inline intT check(nodeT *nd) {
  auto bb = boundingBox<dim>(nd->getItems(), nd->size());
  auto bbn = nd->getBox();

  //check bounding box
  if (bb.first != bbn.first || bb.second != bbn.second) {
    cout << "bounding box error" << endl;
    abort();
  }

  if (!nd->isLeaf()) {
    intT lSize = check<nodeT, dim>(nd->L());
    intT rSize = check<nodeT, dim>(nd->R());

    //check size
    if (lSize + rSize != nd->size()) {
      cout << "size error" << endl;
      abort();
    }

    //check children bounding box
    auto lbb = nd->L()->getBox();
    auto rbb = nd->R()->getBox();
    lbb.first.minCoords(rbb.first);
    lbb.second.maxCoords(rbb.second);
    if (lbb.first != bbn.first || lbb.second != bbn.second) {
      cout << "bounding box inconsistent with children" << endl;
      abort();
    }

  } else {
    // check leaf size
    if (nd->size() != 1) {
      cout << "leaf size error" << endl;
      abort();
    }

    //check leaf bounding box
    if (bbn.first != bbn.second) {
      cout << "wrong bb for leaf" << endl;
      abort();
    }
  }

  return nd->size();
}

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef dynamicKdTree<dim, pointT> dyntreeT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  static const bool serial = false;
  cout << "test kd tree, " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0, t1;

  intT q = 1;
  double range = 5;

  t0.start();
  treeT* tree = new treeT(P, n, !serial, 1);
  cout << ">>>>> range query " << tree->rangeNeighbor(&P[q], range)->size() << endl << endl;
  cout << "static = " << t0.stop() << endl;

  check<nodeT, dim>(tree->rootNode());
  cout << "checked (non-coarsened) tree structure" << endl;

  double init = 0.1;
  intT batches = 10;
  dyntreeT* tree2 = new dyntreeT(P, n*init);
  intT inserted = n*init;
  intT batchSize = (n-inserted)/batches;
  cout << ">>>>> range query " << tree2->rangeNeighbor(&P[q], range)->size() << endl;

  double queryTime = 0;
  t0.start();
  for(int i=0; i<batches; ++i) {
    if(true) {
      tree2->insert(P+inserted, batchSize);
    } else {
      cout << "1: " << tree2->size() << endl << endl;
      tree2->insert(P+inserted, batchSize);
      cout << "2: " << tree2->size() << endl << endl;
      tree2->erase(P+inserted, batchSize);
      cout << "3: " << tree2->size() << endl << endl;
      tree2->insert(P+inserted, batchSize);
      cout << "4: " << tree2->size() << endl << endl;
    }
    inserted += batchSize;
    cout << ">>>>> range query " << tree2->rangeNeighbor(&P[q], range)->size() << endl;
    t1.start();
    par_for(intT i=0; i<inserted/10; ++i) {
      tree2->rangeNeighbor(&P[rand()%inserted], range);
    }
    queryTime += t1.stop();
  }
  cout << "dynamic insert = " << t0.stop() << endl;
  cout << " query time = " << queryTime << endl;

  queryTime = 0;
  t0.start();
  for(int i=0; i<batches; ++i) {
    tree2->erase(P+inserted-batchSize, batchSize);
    inserted -= batchSize;
    cout << ">>>>> range query " << tree2->rangeNeighbor(&P[q], range)->size() << endl;
    t1.start();
    par_for(intT i=0; i<inserted/10; ++i) {
      tree2->rangeNeighbor(&P[rand()%inserted], range);
    }
    queryTime += t1.stop();
  }
  cout << "dynamic erase = " << t0.stop() << endl;
  cout << " query time = " << queryTime << endl;
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
