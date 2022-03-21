#include <gtest/gtest.h>
//#include "common/geometryIO.h"
#include "pargeo/pointIO.h"

//#include <kdtree/binary-heap-layout/bhlkdtree.h>
#include "dynamicKdTree/binary-heap-layout/bhlkdtree.h"

#include "../shared/Shared2DTest.h"
#include "../shared/QueryTest.h"
#include "BHL2DStructureTest.h"

static constexpr int dim = 2;
// <dim, objT, parallel, false>

// not parallel, not coarse
typedef BHL_KdTree<dim, pargeo::point<dim>, false, false> serialSingleTreeT;

INSTANTIATE_TYPED_TEST_SUITE_P(Serial, BHL2DStructureTest, serialSingleTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(Serial_BHL, Shared2DTest, serialSingleTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(Serial_BHL, QueryTest, serialSingleTreeT);

// not parallel, coarse
typedef BHL_KdTree<dim, pargeo::point<dim>, false, true> serialCoarseTreeT;

INSTANTIATE_TYPED_TEST_SUITE_P(SerialCoarse_BHL, Shared2DTest, serialCoarseTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(SerialCoarse_BHL, QueryTest, serialCoarseTreeT);

// parallel, not coarse
typedef BHL_KdTree<dim, pargeo::point<dim>, true, false> parallelSingleTreeT;

INSTANTIATE_TYPED_TEST_SUITE_P(Parallel, BHL2DStructureTest, parallelSingleTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(Parallel_BHL, Shared2DTest, parallelSingleTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(Parallel_BHL, QueryTest, parallelSingleTreeT);

// parallel, coarse
typedef BHL_KdTree<dim, pargeo::point<dim>, true, true> parallelCoarseTreeT;

INSTANTIATE_TYPED_TEST_SUITE_P(ParallelCoarse_BHL, Shared2DTest, parallelCoarseTreeT);
INSTANTIATE_TYPED_TEST_SUITE_P(ParallelCoarse_BHL, QueryTest, parallelCoarseTreeT);
