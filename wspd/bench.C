#include "kdTree.h"
#include "wspd.h"
#include "wspdNormal.h"
#include "pbbs/gettime.h"
using namespace std;

// *************************************************************
//    DRIVER
// *************************************************************

template<int dim>
void bench(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef kdTree<dim, pointT> treeT;
  typedef kdNode<dim, pointT> nodeT;
  static const bool serial = false;
  cout << "test WSPD, " << n << ", dim " << dim << " points" << endl;
  if (n < 2) abort();

  timing t0;
  t0.start();
  treeT* tree = new treeT(P, n, !serial, 1);
  cout << "build-tree-time = " << t0.next() << endl;

  cout << endl << "basic tests" << endl;
  cout << "diag = " << tree->rootNode()->nodeDiag() << endl;
  cout << "node dist = " << tree->rootNode()->L()->nodeDistance(tree->rootNode()->R()) << endl;
  cout << "node far dist = " << tree->rootNode()->L()->nodeFarDistance(tree->rootNode()->R()) << endl;
  cout << "well-separated = " << tree->rootNode()->L()->wellSeparated(tree->rootNode()->R()) << endl;

  auto bcp = tree->rootNode()->L()->compBcp(tree->rootNode()->R());
  cout << "bcp = " << bcp.dist << endl;//todo

  typedef wsp<nodeT> pairT;
  typedef struct nodeT::bcp bcpT;

  t0.stop();t0.start();
  cout << endl << "computes wspd" << endl;
  vector<pairT> out;
  auto wg = wspdNormalSerial<nodeT>(&out);
  wspdSerial<nodeT, wspdNormalSerial<nodeT>>(tree->rootNode(), &wg);
  cout << "#wsp = " << out.size() << endl;
  cout << "serial wspd-time = " << t0.next() << endl;

  t0.stop();t0.start();
  cout << endl << "computes wspd" << endl;
  auto wgpar = wspdNormalParallel<nodeT>(tree->rootNode()->size());
  wspdParallel<nodeT, wspdNormalParallel<nodeT>>(tree->rootNode(), &wgpar);
  auto out2 = wgpar.collect();
  cout << "#wsp = " << out2->size() << endl;
  cout << "parallel wspd-time = " << t0.next() << endl;

  cout << endl << "computes bccp for all pairs" << endl;
  auto bcps = newA(bcpT, out.size());
  par_for (intT i=0; i<out.size(); ++i) {
    auto bcp = out[i].u->compBcp(out[i].v);
    bcps[i] = bcp;
  }
  cout << "bcp-time = " << t0.stop() << endl;

  free(bcps);
}

template void bench<2>(point<2>*, intT);
template void bench<3>(point<3>*, intT);
template void bench<4>(point<4>*, intT);
template void bench<5>(point<5>*, intT);
template void bench<6>(point<6>*, intT);
template void bench<7>(point<7>*, intT);
template void bench<8>(point<8>*, intT);
template void bench<9>(point<9>*, intT);
