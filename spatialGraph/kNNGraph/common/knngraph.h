#include "geometry.h"

template<int dim>
tuple<intT*, intT*> knnGraph(point<dim>* P, intT n, intT k, bool directed);