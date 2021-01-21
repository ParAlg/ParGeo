#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#include "tree.h"
#include "kdt.h"

using namespace std;

template<int dim>
void tree(parlay::sequence<point<dim>> const &Points) {
  size_t n = Points.size();

  bool parallel = true;

  kdNode<dim, point<dim>>* tree = buildKdt<dim, point<dim>>(Points, parallel);

  free(tree);
}

template void tree<2>(parlay::sequence<point<2>> const &);
template void tree<3>(parlay::sequence<point<3>> const &);
