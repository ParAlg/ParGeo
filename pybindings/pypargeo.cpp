//#include "common/get_time.h"
#include "spatialGraph.h"
#include "common/geometryIO.h"
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/*
  Basic wrapper functions
*/

template<class T, class Seq>
py::array_t<T> wrap_array_2d(Seq& result_vec, ssize_t width = 2) {
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = { (ssize_t) result_vec.size(), (ssize_t) width };
  std::vector<ssize_t> strides = { (ssize_t) sizeof(size_t) * (ssize_t) width , (ssize_t) sizeof(size_t) };

  return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				   sizeof(size_t),                          /* size of one scalar        */
				   py::format_descriptor<T>::format(),      /* data type                 */
				   2,                                       /* number of dimensions      */
				   shape,                                   /* shape of the matrix       */
				   strides                                  /* strides for each axis     */
				   ));
}

/*
  IO functions
*/

py::array_t<double> py_loadPoints(std::string& fileName) {
  int dim = readHeader(fileName.c_str());

  if (dim == 2) {
    auto pts = readPointsFromFile<pargeo::point<2>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 2);
  } else if (dim == 3) {
    auto pts = readPointsFromFile<pargeo::point<3>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 3);
  } else if (dim == 4) {
    auto pts = readPointsFromFile<pargeo::point<4>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 4);
  } else if (dim == 5) {
    auto pts = readPointsFromFile<pargeo::point<5>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 5);
  } else if (dim == 6) {
    auto pts = readPointsFromFile<pargeo::point<6>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 6);
  } else if (dim == 7) {
    auto pts = readPointsFromFile<pargeo::point<7>>(fileName.c_str());
    return wrap_array_2d<double>(pts, 7);
  } else throw std::runtime_error("dimensionlity not yet supported");
}

/*
  Spatial graph generation wrappers
*/

py::array_t<size_t> py_delaunayGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
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
    return wrap_array_2d<size_t>(result_vec);

  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<size_t> py_gabrielGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array) {
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
    return wrap_array_2d<size_t>(result_vec);

  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<size_t> py_skeleton(py::array_t<double, py::array::c_style | py::array::forcecast> array, double beta) {
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
    return wrap_array_2d<size_t>(result_vec);

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

  if (dim == 2) {

    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<2>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else if (dim == 3) {

    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<3>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else if (dim == 4) {

    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<4>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else if (dim == 5) {

    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<5>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else if (dim == 6) {

    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<6>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else if (dim == 7) {

    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<7>(array_vec, k);
    return wrap_array_2d<size_t>(result_vec);

  } else
    throw std::runtime_error("Only supports 2-7d points.");
}

PYBIND11_MODULE(pypargeo, m)
{
  m.doc() = "Pargeo: library for parallel computational geometry.";

  m.def("loadPoints", &py_loadPoints, "Load points based on file path.");

  m.def("KnnGraph", &py_knnGraph, "Input: 2d-numpy array containing n points in R^2; parameter k. Output: array E of size n*k, where E[i] and E[i+1] forall i are point indices represents a directed edge of the directed-knn graph.");

  m.def("DelaunayGraph", &py_delaunayGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of edges, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the delaunay graph.");

  m.def("GabrielGraph", &py_gabrielGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the Gabriel graph.");

  m.def("BetaSkeleton", &py_skeleton, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the beta skeleton graph.");
}
