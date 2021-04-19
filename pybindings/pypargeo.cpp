//#include "common/get_time.h"
#include "spatialGraph.h"
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/*
  Basic wrapper functions
*/

template<class T, class Seq>
py::array_t<T> wrap_array_2d(Seq& result_vec) {
  ssize_t ndim = 2;
  std::vector<ssize_t> shape = { (ssize_t) result_vec.size(), 2 };
  std::vector<ssize_t> strides = { sizeof(size_t) * 2 , sizeof(size_t) };

  return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				   sizeof(size_t),                          /* size of one scalar        */
				   py::format_descriptor<T>::format(),      /* data type                 */
				   2,                                       /* number of dimensions      */
				   shape,                                   /* shape of the matrix       */
				   strides                                  /* strides for each axis     */
				   ));
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

  } else
    throw std::runtime_error("Only supports 2d points.");
}

PYBIND11_MODULE(pypargeo, m)
{
  m.doc() = "Pargeo: library for parallel computational geometry.";

  m.def("KnnGraph", &py_knnGraph, "Input: 2d-numpy array containing n points in R^2; parameter k. Output: array E of size n*k, where E[i] and E[i+1] forall i are point indices represents a directed edge of the directed-knn graph.");

  m.def("DelaunayGraph", &py_delaunayGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of edges, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the delaunay graph.");

  m.def("GabrielGraph", &py_gabrielGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the Gabriel graph.");

  m.def("BetaSkeleton", &py_skeleton, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the beta skeleton graph.");
}
