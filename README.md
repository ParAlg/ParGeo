# **  Forked ParGeo: A Library for Parallel Computational Geometry**

ParGeo is a library containing a collection of implementations of parallel algorithms in computational geometry in low-dimensional space.

ParGeo contains efficient multicore implementations of kd-trees. 
The code supports kd-tree based spatial search, including
K-NN and range search.
The code can also compute the bichromatic closest pair by traversing two kd-trees.
The code is optimized for fast kd-tree construction by performing the split inparallel,
and the queries themselves are data-parallel.
We also include a parallel batch-dynamic kd-tree that supports batch insertions and
deletions.
Our kd-tree can be used to generate a well-separated pair
decomposition,
which can be used to compute the Euclidean minimum spanning tree (EMST) and spanners.

In addition, ParGeo contains parallel implementations for classic algorithms in computational
geometry, including Morton sorting, closest pair, convex hull, and smallest enclosing ball.
ParGeo also contains a collection of geometric graph generators
for point data sets.
It includes routines for constructing K-NN graphs using kd-trees,
and also supports common spatial network graphs, including the
Delaunay graph and the beta-skeleton graph.

## **Content**

- [**ParGeo: A Library for Parallel Computational Geometry**](#pargeo-a-library-for-parallel-computational-geometry)
  - [**Content**](#content)
  - [**Code Organization**](#code-organization)
  - [**Tutorial**](#tutorial)
    - [**Getting Started**](#getting-started)
    - [**Multi-core Parallelism**](#multi-core-parallelism)
    - [**Benchmarking**](#benchmarking)
    - [**Tests**](#tests)
  - [**Modules**](#modules)
    - [**Data Point Representation and IO**](#data-point-representation-and-io)
    - [**Kd-Tree**](#kd-tree)
      - [**K-NN Search**](#k-nn-search)
      - [**Range Search**](#range-search)
      - [**Well-Separated Pair Decomposition (WSPD)**](#well-separated-pair-decomposition-wspd)
      - [**Bichromatic Closest Pair**](#bichromatic-closest-pair)
      - [**Dynamic Kd-Tree**](#dynamic-kd-tree)
    - [**Spatial Graph Generator**](#spatial-graph-generator)
      - [**K-NN Graph**](#k-nn-graph)
      - [**Delaunay Graph**](#delaunay-graph)
      - [**Beta-Skeleton and Gabriel Graph**](#beta-skeleton-and-gabriel-graph)
      - [**T-Spanner**](#t-spanner)
    - [**Spatial Clustering**](#spatial-clustering)
    - [**Convex Hull**](#convex-hull)
    - [**Smallest Enclosing Ball**](#smallest-enclosing-ball)
    - [**Closest Pair**](#closest-pair)
    - [**Morton Sort**](#morton-sort)
    - [**Data Generator**](#data-generator)
      - [**Uniform Data Generator**](#uniform-data-generator)
      - [**Clustered Data Generator**](#clustered-data-generator)
  - [**References**](#references)

## **Code Organization**

- **benchmark/** Code for benchmarking and comparing implementations using [Google Benchmark](https://github.com/google/benchmark).
- **example/** Examples and executable files for the modules in ParGeo,
- **external/** External dependencies (git submodules).
- **include/** Header files.
- **src/** Source files.
- **test/** Tests using [Google Test](https://github.com/google/googletest).

While some of the implementations are header-only, most of the modules have a folder in both the `include/` and `src/` folder, containing header and source respectively.

## **Tutorial**

### **Getting Started**

ParGeo uses the [parlaylib](https://github.com/cmuparlay/parlaylib) as a submodule for multi-core parallel programming. Before building the project, initialize the submodules:

```
git submodule init
git submodule update
```

ParGeo uses CMake. To build the entire project:

```
mkdir build
cd build
cmake ..
make -j
```

After the build completes, you will be able to find the executables in the nested directories of `build/`. Running `./program` in the terminal will return short usage instructions. A typical run will involve `./program <input-arguments> <data-path>`.

You don't have to build the entire project, which can be very slow. For example, to only build the spatial graph generator:

```
mkdir build
cd build/example
cmake ..
make -j graphGenerator
```

To generate a 1-NN graph:

```
./graphGenerator -algo 1 -param 1 -o edges.txt dataset.txt
```

ParGeo automatically parses csv-like point data files with or without header. The delimiters (space, comma etc) have to be one character each in length. [Here](https://github.com/ParAlg/ParGeo/blob/main/test/datasets/2d_100.txt) is an example of a 2-dimensional data set.

[back to top](#content)

### **Multi-core Parallelism**

### **Benchmarking**

The `benchmark/` folder contains benchmarking code for the modules in ParGeo using [Google Benchmark](https://github.com/google/benchmark). The benchmark generates data set at runtime using ParGeo's data generator (`include/dataset/`).

[back to top](#content)

### **Tests**

The `test/` folder contains test code for the modules in ParGeo using [Google Test](https://github.com/google/googletest). To run all the tests, first compile the tests, and then run `ctest` from the `build/test` directory.

[back to top](#content)

## **Modules**

### **Data Point Representation and IO**

```
#include "pargeo/point.h"
```

`include/pargeo/point.h` contains the header-only implementations of data point classes used throughout ParGeo. The dimensionality of the data set is passed in as a template parameter, therefore the size of each instantiation of the point is the dimensionality times the size of each data field. For example, as defined towards the bottom of `point.h`, `point<2>` has the size of two doubles, while `fpoint<3>` has the size of three floats. The functions that use the point classes will need to either define their dimensionality explicitly, or be templated. Most implementations in ParGeo have dimensionality defined for 2 to 7, and they can be extended to other dimensionality by adding additional definitions.

The point classes are augmented by a number of subroutines to access, modify the data fields, and perform basic arithmetic and geometric computations. The IO of the points are implemented in `include/pargeo/pointIO.h`. Examples of data IO can be found in most of the examples in `example/`.

### **Kd-Tree**

```
#include "kdTree/kdTree.h"
```

The kd-tree is a spatial indexing data structure that recursively divide the data space into two halves using axis-aligned hyperplanes, forming a binary tree. ParGeo contains a header-only implementation of a parallel kd-tree, and a number of related algorithms detailed below.

[back to top](#content)

#### **K-NN Search**

```
#include "kdTree/kdTree.h"
```

K-NN search is the problem of finding the k-nearest points in a point set given a query point that may or may not be in the set. ParGeo's header-only implementation in `include/kdTree/knnImpl.h` uses the kd-tree implementation in `include/kdTree/kdTree.h` to find the k-nearest points under the Euclidean distance metric. ParGeo uses a classic depth-first search algorithm while keeping the top-k points during the traversal. Examples of using the K-NN search can be found in `example/kNearestNeighbor.cpp`.

[back to top](#content)

#### **Range Search**

```
#include "kdTree/kdTree.h"
```

Range search is the problem of finding the points that intersects a query object. ParGeo supports two commonly used query objects, hypersphere and hyperrectangle. ParGeo's header-only implementation of the range search in `include/kdTree/rangeSearchImpl.h` is built on the kd-tree implementation in `include/kdTree/kdTree.h`.
Examples of using the range search can be found in `example/rangeSearch.cpp`.

[back to top](#content)

#### **Well-Separated Pair Decomposition (WSPD)**

```
#include "kdTree/kdTree.h"
```

The [well-separated pair decomposition](https://en.wikipedia.org/wiki/Well-separated_pair_decomposition) (WSPD) of a point set is a sequence of pairs of sets, such that each pair is well-separated. For each two distinct points, there exists one and exactly one pair that separates them. ParGeo contains a header-only implementation of the WSPD in `include/kdTree/wspdImpl.h`, based on the fair-split kd-tree in `include/kdTree/kdTree.h`. Examples of computing the WSPD can be found in `example/wellSeparatedPair.cpp`.

[back to top](#content)

#### **Bichromatic Closest Pair**

```
#include "kdTree/kdTree.h"
```

Bichromatic closest pair (BCP) is the problem of finding the closest pair of distinctly colored points given two sets of points of two different colors. The naive algorithm for the BCP compares all pairs of points between the sets. ParGeo uses a more efficient implementation by contructing a kd-tree on each set, and traverse the two trees while maintaining the closest pair, and pruning nodes further than the closest pair.
ParGeo contains a header-only implementation of the algorithm in `include/kdTree/bccpImpl.h`, based on the tree in `include/kdTree/kdTree.h`. Examples of computing the BCP can be found in `example/bichromaticClosestPair.cpp`.

[back to top](#content)

#### **Dynamic Kd-Tree**

Currently ParGeo contains a simple dynamic kd-tree that supports batch-insertion and batch-deletion of items. It also supports k-nearest neighbor search. The implementation is in `include/kdTree/dynKdTree.h`, and the usage examples is in `dynamicKdTreeKNN.cpp`.

[back to top](#content)

### **Spatial Graph Generator**

```
#include "spatialGraph/spatialGraph.h"
```

ParGeo contains efficient parallel implementations of a number of spatial graph generators. A binary executable, and examples for all the generators can be found in `example/graphGenerator.cpp`. The generators are described in greater details in [[2](https://people.csail.mit.edu/jshun/geograph.pdf)], and we provide some brief descritions below.

#### **K-NN Graph**

ParGeo generates the K-NN graph on a point data set, and outputs edges of a directed graph. The K-NN graph is generated by computing the k-nearest neighbors of all points using the kd-tree. Given a data set of size n, for example, if k = 1, there will be n edges.

#### **Delaunay Graph**

The Delaunay graph is a graph defined by the edges of the [delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation). ParGeo generates the undirected Delaunay graph of 2d point data sets. It makes use of an implementation of Delaunay triangulation of the [pbbsbench](https://github.com/cmuparlay/pbbsbench).

#### **Beta-Skeleton and Gabriel Graph**

The [Beta skeleton](https://en.wikipedia.org/wiki/Beta_skeleton) is an undirected graph defined from a set of points on the Euclidean plane. ParGeo implements the lune-based beta-skeleton as defined. The Gabriel graph is a special case of the beta skeleton where beta = 1. ParGeo implements the beta skeleton by removing edges from a Delaunay graph. Specifically, an edge is removed if the lune is empty. ParGeo determines the emptiness of the lune using range searches supported by a kd-tree.

#### **T-Spanner**

The [t-spanner](https://en.wikipedia.org/wiki/Geometric_spanner) graph is an undirected graph where there is a t-path between any pair of vertices for the parameter t. A t-path is defined as a path through the graph with weight at most t-times the spatial distance between its endpoints. ParGeo contains an Euclidean t-spanner generator based on the well-separated pair decomposition (WSPD). Specifically, an edge is added between an arbitrary pair of points between every well-separated pair.

[back to top](#content)

### **Spatial Clustering**

TBA

[back to top](#content)

### **Convex Hull**

```
#include "convexHull2d/divideConquer/hull.h"
#include "convexHull3d/divideConquer/hull.h"
```

The convex hull of a point set P is the minimal [convex set](https://en.wikipedia.org/wiki/Convex_hull) containing P. ParGeo contains a number of convex hull implementations , for both 2 and 3 dimensional data sets. For both 2 and 3 dimensions, we use a similar parallelization strategy. Specifically, we block the data into number of blocks proportional to the number of processors. The convex hull for each block is computed, and then a final convex hull is computed on the vertices of the sub-convex hulls. We also explore and evaluate other convex hull algorithms, defailed in [[1](https://www.cs.ucr.edu/~ygu/papers/PPoPP22/ParGeo.pdf)].

The implementation for 2 and 3 dimension can be found in `include/convexHull2d/divideConquer/hull.h` and `include/convexHull3d/divideConquer/hull.h` respectively. For both implementations, the algorithm that computes the each sub-convex hull sequentially is the [quickhull](https://en.wikipedia.org/wiki/Quickhull) algorithm. Examples of running convex hull algorithms can be found in `example/convexHull.cpp`.

[back to top](#content)

### **Smallest Enclosing Ball**

```
#include "enclosingBall/sampling/seb.h"
```

The [smallest enclosing ball](https://en.wikipedia.org/wiki/Smallest-circle_problem) (SEB) is the problem of finding the smallest enclosing bounding sphere of a given set of points. ParGeo contains a number of SEB implementations, among which `include/enclosingBall/hull.h` is the fastest. The algorithm is based on orthant-scan, but with a sampling technique described in [[1](https://www.cs.ucr.edu/~ygu/papers/PPoPP22/ParGeo.pdf)]. Examples of running convex hull algorithms can be found in `example/smallestEnclosingBall.cpp`.

[back to top](#content)

### **Closest Pair**

```
#include "closestPair/closestPair.h"
```

The [closest pair](https://en.wikipedia.org/wiki/Closest_pair_of_points_problem) is the problem of finding a pair of points with the smallest distance between them. ParGeo implements a parallel divide-and-conquer algorithm proposed by Blelloch et. al. An usage example can be found in `example/closestPair.cpp`.

[back to top](#content)

### **Morton Sort**

```
#include "mortonSort/mortonSort.h"
```

ParGeo performs [Morton sort](https://en.wikipedia.org/wiki/Z-order_curve) of 2 and 3 dimensional point sets. Usage examples can be found in `example/mortonSort.cpp`.

[back to top](#content)

### **Data Generator**

ParGeo contains various synthetic data generators that produce synthetic data of various density.

#### **Uniform Data Generator**

```
#include "dataset/uniform.h"
```

ParGeo generates uniformly distributed points constrained by either a polytope, either a hypersphere or a hypercube. User also has the optional of specifying whether the points should be distributed within the polytope, or on the polytope's faces. A uniform data generator can be found at `example/uniformDataGenerator.cpp`.

#### **Clustered Data Generator**

```
#include "dataset/seedSpreader.h"
```

ParGeo implements the random seed spreader proposed by [Gan and Tao](https://dl.acm.org/doi/10.1145/3083897). The seed spreader generates random clustered point data via a random walk. A seed spreader executable can be found at `example/seedSpreaderGenerator.cpp`. User can choose whether to generate clusters with similar density or variable densities.

[back to top](#content)

## **References**

[1] [Yiqiu Wang, Shangdi Yu, Laxman Dhulipala, Yan Gu, Julian Shun. POSTER: ParGeo: A Library for Parallel Computational Geometry. PPoPP 2022.](https://www.cs.ucr.edu/~ygu/papers/PPoPP22/ParGeo.pdf)

[2] [Yiqiu Wang, Shangdi Yu, Laxman Dhulipala, Yan Gu, Julian Shun. GeoGraph: A Framework for Graph Processing on Geometric Data. SIGOPS Operating Systems Review 55.](https://people.csail.mit.edu/jshun/geograph.pdf)
