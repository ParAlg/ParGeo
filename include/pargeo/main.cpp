#include "dynKdTree.h"

template<int dim>
class point: public dynKdTree::coordinate<dim> {

public:

  point(double* _data): dynKdTree::coordinate<dim>(_data) {

  }

};


int main() {

  using namespace dynKdTree;

  return 0;
}
