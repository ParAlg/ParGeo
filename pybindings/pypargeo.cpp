#ifndef SILENT
#define VERBOSE
#endif

#ifdef VERBOSE
#include "common/get_time.h"
#endif

#include "spatialGraph.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<size_t> py_delaunay(py::array_t<double, py::array::c_style | py::array::forcecast> array)
{

  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];

  size_t n = array.size() / dim;

#ifdef VERBOSE
  timer t;
#endif

  if (dim == 2) {
#ifdef VERBOSE
    t.start();
#endif
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
#ifdef VERBOSE
    std::cout << "::copy-in-time = " << t.get_next() << std::endl;
#endif

    parlay::sequence<pargeo::edge> result_vec = pargeo::delaunayGraph<2>(array_vec);

#ifdef VERBOSE
    std::cout << "::compute-time = " << t.get_next() << std::endl;
#endif

    ssize_t ndim = 2;
    std::vector<ssize_t> shape = { (ssize_t) result_vec.size(), 2 };
    std::vector<ssize_t> strides = { sizeof(size_t) * 2 , sizeof(size_t) };

    // return a 2d numpy array
    return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				     sizeof(size_t),                          /* size of one scalar        */
				     py::format_descriptor<size_t>::format(), /* data type                 */
				     2,                                       /* number of dimensions      */
				     shape,                                   /* shape of the matrix       */
				     strides                                  /* strides for each axis     */
				     ));

  } else {
    throw std::runtime_error("Only supports 2d points.");
  }

}

py::array_t<size_t> py_gabriel(py::array_t<double, py::array::c_style | py::array::forcecast> array)
{

  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];

  size_t n = array.size() / dim;

#ifdef VERBOSE
  timer t;
#endif

  if (dim == 2) {
#ifdef VERBOSE
    t.start();
#endif
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
#ifdef VERBOSE
    std::cout << "::copy-in-time = " << t.get_next() << std::endl;
#endif

    parlay::sequence<pargeo::edge> result_vec = pargeo::gabrielGraph<2>(array_vec);

#ifdef VERBOSE
    std::cout << "::compute-time = " << t.get_next() << std::endl;
#endif

    ssize_t ndim = 2;
    std::vector<ssize_t> shape = { (ssize_t)result_vec.size(), 2 };
    std::vector<ssize_t> strides = { sizeof(size_t) * 2 , sizeof(size_t) };

    // return a 2d numpy array
    return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				     sizeof(size_t),                          /* size of one scalar        */
				     py::format_descriptor<size_t>::format(), /* data type                 */
				     2,                                       /* number of dimensions      */
				     shape,                                   /* shape of the matrix       */
				     strides                                  /* strides for each axis     */
				     ));

  } else {
    throw std::runtime_error("Only supports 2d points.");
  }

}

py::array_t<size_t> py_skeleton(py::array_t<double, py::array::c_style | py::array::forcecast> array, double beta)
{

  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];

  size_t n = array.size() / dim;

#ifdef VERBOSE
  timer t;
#endif

  if (dim == 2) {
#ifdef VERBOSE
    t.start();
#endif
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
#ifdef VERBOSE
    std::cout << "::copy-in-time = " << t.get_next() << std::endl;
#endif

    parlay::sequence<pargeo::edge> result_vec = pargeo::betaSkeleton<2>(array_vec, beta);

#ifdef VERBOSE
    std::cout << "::compute-time = " << t.get_next() << std::endl;
#endif

    ssize_t ndim = 2;
    std::vector<ssize_t> shape = { (ssize_t)result_vec.size(), 2 };
    std::vector<ssize_t> strides = { sizeof(size_t) * 2 , sizeof(size_t) };

    // return a 2d numpy array
    return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				     sizeof(size_t),                          /* size of one scalar        */
				     py::format_descriptor<size_t>::format(), /* data type                 */
				     2,                                       /* number of dimensions      */
				     shape,                                   /* shape of the matrix       */
				     strides                                  /* strides for each axis     */
				     ));

  } else {
    throw std::runtime_error("Only supports 2d points.");
  }

}

py::array_t<size_t> py_knngraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k)
{

  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if (sizeof(pargeo::point<2>) != 16)
    throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

  int dim = array.shape()[1];

  size_t n = array.size() / dim;

#ifdef VERBOSE
  timer t;
#endif

  if (dim == 2) {
#ifdef VERBOSE
    t.start();
#endif
    parlay::sequence<pargeo::point<2>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
#ifdef VERBOSE
    std::cout << "::copy-in-time = " << t.get_next() << std::endl;
#endif

    parlay::sequence<pargeo::edge> result_vec = pargeo::knnGraph<2>(array_vec, k);

#ifdef VERBOSE
    std::cout << "::compute-time = " << t.get_next() << std::endl;
#endif

    ssize_t ndim = 2;
    std::vector<ssize_t> shape = { (ssize_t)result_vec.size(), 2 };
    std::vector<ssize_t> strides = { sizeof(size_t) * 2 , sizeof(size_t) };

    // return a 2d numpy array
    return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				     sizeof(size_t),                          /* size of one scalar        */
				     py::format_descriptor<size_t>::format(), /* data type                 */
				     2,                                       /* number of dimensions      */
				     shape,                                   /* shape of the matrix       */
				     strides                                  /* strides for each axis     */
				     ));
  } else {
    throw std::runtime_error("Only supports 2d points.");
  }
}

PYBIND11_MODULE(pypargeo, m)
{
  m.doc() = "Generates spatial graphs for a point data set in R^2.";

  m.def("knn_graph", &py_knngraph, "Input: 2d-numpy array containing n points in R^2; parameter k. Output: array E of size n*k, where E[i] and E[i+1] forall i are point indices represents a directed edge of the directed-knn graph.");

  m.def("delaunay_graph", &py_delaunay, "Input: 2d-numpy array containing n points in R^2. Output: array E of edges, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the delaunay graph.");

  m.def("gabriel_graph", &py_gabriel, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the Gabriel graph.");

  m.def("beta_skeleton", &py_skeleton, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the beta skeleton graph.");
}
