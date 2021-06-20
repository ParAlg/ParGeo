#include "spatialGraph/spatialGraph.h"
#include "euclideanMst/euclideanMst.h"
#include "pargeo/pointIO.h"
#include "pargeo/graphIO.h"
#include <string>
#include <limits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pargeo::graphIO;
using namespace pargeo::pointIO;

/*
  Helper functions
*/

template<class T, class Seq>
py::array_t<T> wrapArray2d(Seq& result_vec, ssize_t cols) {
  ssize_t rows = (ssize_t) result_vec.size() *
    (ssize_t) sizeof(typename Seq::value_type) /
    (ssize_t) sizeof(T);
  rows /= cols;
  std::vector<ssize_t> shape = { rows, cols };
  std::vector<ssize_t> strides = { (ssize_t) sizeof(T) * cols, (ssize_t) sizeof(T) };

  return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				   sizeof(T),                               /* size of one scalar        */
				   py::format_descriptor<T>::format(),      /* data type                 */
				   2,                                       /* number of dimensions      */
				   shape,                                   /* shape of the matrix       */
				   strides                                  /* strides for each axis     */
				   ));
}

template<class T, class Edg>
py::array_t<T> castEdgeList(Edg& E) {
  sequence<T> result_vec(E.size()*2);
  parlay::parallel_for(0, E.size(),
	       [&](size_t i){
		 result_vec[i*2] = (T) E[i].u;
		 result_vec[i*2+1] = (T) E[i].v;
	       });
  return wrapArray2d<T>(result_vec, 2);
}

// Computes edge weights;
// removes edges with weight >epsilon
template<class T, class Edg, class Pts>
py::array_t<T> castWghEdgeList(Edg& E,
			       Pts& P,
			       T eps=-1) {
  struct wghEdge {T u, v, weight;};
  sequence<wghEdge> result(E.size());
  parlay::parallel_for(0, E.size(),
	       [&](size_t i){
		 result[i].u = (T) E[i].u;
		 result[i].v = (T) E[i].v;
		 result[i].weight = (T) P[E[i].u].dist(P[E[i].v]);
	       });
  if (eps > 0) {
    auto filtered =
      parlay::filter(result, [&](auto e){ return e.weight <= eps; });
    return wrapArray2d<T>(filtered, 3);
  } else {
    return wrapArray2d<T>(result, 3);
  }
}

template<class T, class Edg>
py::array_t<T> castWghEdgeList(Edg& E,
			       T eps=-1) {
  struct wghEdge {T u, v, weight;};
  sequence<wghEdge> result(E.size());
  parlay::parallel_for(0, E.size(),
	       [&](size_t i){
		 result[i].u = (T) E[i].u;
		 result[i].v = (T) E[i].v;
		 result[i].weight = (T) E[i].weight;
	       });
  if (eps > 0) {
    auto filtered =
      parlay::filter(result, [&](auto e){ return e.weight <= eps; });
    return wrapArray2d<T>(filtered, 3);
  } else {
    return wrapArray2d<T>(result, 3);
  }
}

/*
  IO functions
*/

/* to unweighted SNAP format on disk */
void py_writeEdges(std::string& fileName, py::array_t<size_t, py::array::c_style | py::array::forcecast> array) {
  // todo check the type of array

  int dim = array.shape()[1]; // number of entries per line
  size_t n = array.size() / dim;

  if (dim == 2) {
    struct edge {size_t u, v;};
    parlay::sequence<edge> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(size_t));
    writeEdgeSeqToFile<edge>(array_vec, fileName.c_str());
  } else
    throw std::runtime_error("Invalid unweighted edge array, can't save to disk.");
}

/* to SNAP format on disk */
void py_writeWghEdges(std::string& fileName, py::array_t<float, py::array::c_style | py::array::forcecast> array) {
  // todo check the type of array

  int dim = array.shape()[1]; // number of entries per line
  size_t n = array.size() / dim;

  if (dim == 3) {
    struct edge {float u, v, weight;};
    parlay::sequence<edge> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(float));
    writeWghEdgeSeqToFile<edge>(array_vec, fileName.c_str());
  } else
    throw std::runtime_error("Invalid weighted python edge array, can't save to disk.");
}

py::array_t<double> py_loadPoints(std::string& fileName) {
  int dim = readHeader(fileName.c_str());

  if (dim == 2) {
    auto pts = readPointsFromFile<pargeo::point<2>>(fileName.c_str());
    return wrapArray2d<double>(pts, 2);
  } else if (dim == 3) {
    auto pts = readPointsFromFile<pargeo::point<3>>(fileName.c_str());
    return wrapArray2d<double>(pts, 3);
  } else if (dim == 4) {
    auto pts = readPointsFromFile<pargeo::point<4>>(fileName.c_str());
    return wrapArray2d<double>(pts, 4);
  } else if (dim == 5) {
    auto pts = readPointsFromFile<pargeo::point<5>>(fileName.c_str());
    return wrapArray2d<double>(pts, 5);
  } else if (dim == 6) {
    auto pts = readPointsFromFile<pargeo::point<6>>(fileName.c_str());
    return wrapArray2d<double>(pts, 6);
  } else if (dim == 7) {
    auto pts = readPointsFromFile<pargeo::point<7>>(fileName.c_str());
    return wrapArray2d<double>(pts, 7);
  } else throw std::runtime_error("dimensionlity not yet supported");
}

/*
  Spatial graph generation wrappers
*/

py::array_t<unsigned int> py_delaunayGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::delaunayGraph<2>(array_vec);
    return castEdgeList<unsigned int>(result_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<float> py_wghDelaunayGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::delaunayGraph<2>(array_vec);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<unsigned int> py_gabrielGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::gabrielGraph<2>(array_vec);
    return castEdgeList<unsigned int>(result_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<float> py_wghGabrielGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::gabrielGraph<2>(array_vec);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<unsigned int> py_skeleton(py::array_t<double, py::array::c_style | py::array::forcecast> array, double beta) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::betaSkeleton<2>(array_vec, beta);
    return castEdgeList<unsigned int>(result_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<float> py_wghSkeleton(py::array_t<double, py::array::c_style | py::array::forcecast> array, double beta) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::betaSkeleton<2>(array_vec, beta);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<size_t> py_knnGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;
  parlay::sequence<pargeo::dirEdge> result_vec;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<2>(array_vec, k);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<3>(array_vec, k);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<4>(array_vec, k);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<5>(array_vec, k);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<6>(array_vec, k);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::knnGraph<7>(array_vec, k);
  } else
    throw std::runtime_error("Only dimensions 2-7 are supported");

  return castEdgeList<unsigned int>(result_vec);
}


py::array_t<float> py_wghKnnGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k, double eps = -1) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;
  parlay::sequence<pargeo::dirEdge> result_vec;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<2>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<3>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<4>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<5>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<6>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::knnGraph<7>(array_vec, k);
    return castWghEdgeList<float>(result_vec, array_vec, eps);
  } else
    throw std::runtime_error("Only dimensions 2-7 are supported");
}

py::array_t<size_t> py_spanner(py::array_t<double, py::array::c_style | py::array::forcecast> array, double t) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;
  parlay::sequence<pargeo::edge> result_vec;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<2>(array_vec, t);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<3>(array_vec, t);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<4>(array_vec, t);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<5>(array_vec, t);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<6>(array_vec, t);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::spanner<7>(array_vec, t);
  } else
    throw std::runtime_error("Only dimensions 2-7 are supported");

  return castEdgeList<unsigned int>(result_vec);
}

py::array_t<float> py_wghSpanner(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t t) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;
  parlay::sequence<pargeo::edge> result_vec;

  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<2>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<3>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<4>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<5>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<6>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    auto result_vec = pargeo::spanner<7>(array_vec, t);
    return castWghEdgeList<float>(result_vec, array_vec);
  } else
    throw std::runtime_error("Only dimensions 2-7 are supported");
}

py::array_t<float> py_wghEmst(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];
  size_t n = array.size() / dim;

  parlay::sequence<pargeo::wghEdge> result_vec;
  if (dim == 2) {
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<2>(array_vec);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<3>(array_vec);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<4>(array_vec);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<5>(array_vec);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<6>(array_vec);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    result_vec = pargeo::euclideanMst<7>(array_vec);
  } else
    throw std::runtime_error("Only supports 2-7 dim points.");

    return castWghEdgeList<float>(result_vec);
}

PYBIND11_MODULE(pypargeo, m)
{
  m.doc() = "Pargeo: a library for parallel computational geometry.";

  m.def("writeEdges",
	&py_writeEdges,
	"Save unweighted edge list in SNAP format to the file path.",
	py::arg("fileName"), py::arg("array"));

  m.def("writeWghEdges",
	&py_writeWghEdges,
	"Save weighted edge list in SNAP format to the file path.",
	py::arg("fileName"), py::arg("array"));

  m.def("loadPoints",
	&py_loadPoints,
	"Load points based on file path.");

  m.def("KnnGraph",
	&py_knnGraph,
	"Unweighted knn-graph of points in R^{2-7}.",
	py::arg("array"), py::arg("k"));

  m.def("WghKnnGraph",
	&py_wghKnnGraph,
	"Weighted knn-graph of points in R^{2-7}.",
	py::arg("array"), py::arg("k"), py::arg("eps")=-1);

  m.def("DelaunayGraph",
	&py_delaunayGraph,
	"Unweighted planar delaunay graph.",
	py::arg("array"));

  m.def("WghDelaunayGraph",
	&py_wghDelaunayGraph,
	"Weighted planar delaunay graph.",
	py::arg("array"));

  m.def("GabrielGraph",
	&py_gabrielGraph,
	"Unweighted planar Gabriel graph.",
	py::arg("array"));

  m.def("WghGabrielGraph",
	&py_wghGabrielGraph,
	"Weighted planar Gabriel graph.",
	py::arg("array"));

  m.def("BetaSkeleton",
	&py_skeleton,
	"Unweighted planar Beta skeleton.",
	py::arg("array"), py::arg("beta"));

  m.def("WghBetaSkeleton",
	&py_wghSkeleton,
	"Weighted planar Beta skeleton.",
	py::arg("array"), py::arg("beta"));

  m.def("Spanner",
	&py_spanner,
	"Unweighted spanner based on the WSPD.",
	py::arg("array"), py::arg("t"));

  m.def("WghSpanner",
	&py_spanner,
	"Weighted spanner based on the WSPD.",
	py::arg("array"), py::arg("t"));

  m.def("WghEuclideanMst",
	&py_wghEmst,
	"Weighted Euclidean MST.",
	py::arg("array"));
}
