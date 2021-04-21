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
  std::vector<ssize_t> shape = { (ssize_t) result_vec.size(), (ssize_t) width };
  std::vector<ssize_t> strides = { (ssize_t) sizeof(T) * (ssize_t) width, (ssize_t) sizeof(T) };

  return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				   sizeof(T),                               /* size of one scalar        */
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

static constexpr double weightScaleFactor = 1;

py::array_t<unsigned int> py_delaunayGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, bool weighted=false) {
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

    if (!weighted) {
      return wrap_array_2d<unsigned int>(result_vec, 2);
    } else {
      struct we { unsigned int u, v, weight; };
      parlay::sequence<we> E(result_vec.size());
      parlay::parallel_for(0, result_vec.size(),
			   [&](size_t i){
			     E[i].u = result_vec[i].u;
			     E[i].v = result_vec[i].v;
			     E[i].weight = (unsigned int) weightScaleFactor*array_vec[E[i].u].dist(array_vec[E[i].v]);
			   });
      return wrap_array_2d<unsigned int>(E, 3);
    }

  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<size_t> py_gabrielGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, bool weighted=false) {
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

    if (!weighted) {
      return wrap_array_2d<unsigned int>(result_vec);
    } else {
      struct we { unsigned int u, v, weight; };
      parlay::sequence<we> E(result_vec.size());
      parlay::parallel_for(0, result_vec.size(),
			   [&](size_t i){
			     E[i].u = result_vec[i].u;
			     E[i].v = result_vec[i].v;
			     E[i].weight = (unsigned int) weightScaleFactor*array_vec[E[i].u].dist(array_vec[E[i].v]);
			   });
      return wrap_array_2d<unsigned int>(E, 3);
    }
  } else
    throw std::runtime_error("Only supports 2d points.");
}

py::array_t<size_t> py_skeleton(py::array_t<double, py::array::c_style | py::array::forcecast> array, double beta, bool weighted) {
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

    if (!weighted) {
      return wrap_array_2d<unsigned int>(result_vec);
    } else {
      struct we { unsigned int u, v, weight; };
      parlay::sequence<we> E(result_vec.size());
      parlay::parallel_for(0, result_vec.size(),
			   [&](size_t i){
			     E[i].u = result_vec[i].u;
			     E[i].v = result_vec[i].v;
			     E[i].weight = (unsigned int) weightScaleFactor*array_vec[E[i].u].dist(array_vec[E[i].v]);
			   });
      return wrap_array_2d<unsigned int>(E, 3);
    }
  } else
    throw std::runtime_error("Only supports 2d points.");
}

template <int dim>
inline py::array_t<size_t> py_knnHelper(parlay::sequence<pargeo::point<dim>> &array_vec, size_t k, bool weighted) {
  parlay::sequence<pargeo::dirEdge> result_vec = pargeo::knnGraph<dim>(array_vec, k);

  if (!weighted) {
    return wrap_array_2d<unsigned int>(result_vec, 2);
  } else {
    struct we { unsigned int u, v, weight; };
    parlay::sequence<we> E(result_vec.size());
    parlay::parallel_for(0, result_vec.size(),
			 [&](size_t i){
			   E[i].u = result_vec[i].u;
			   E[i].v = result_vec[i].v;
			   E[i].weight = (unsigned int) weightScaleFactor*
			     array_vec[E[i].u].dist(array_vec[E[i].v]);
			 });
    return wrap_array_2d<unsigned int>(E, 3);
  }
}

py::array_t<size_t> py_knnGraph(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k, bool weighted=false) {
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
    return py_knnHelper<2>(array_vec, k, weighted);
  } else if (dim == 3) {
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    return py_knnHelper<3>(array_vec, k, weighted);
  } else if (dim == 4) {
    parlay::sequence<pargeo::point<4>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    return py_knnHelper<4>(array_vec, k, weighted);
  } else if (dim == 5) {
    parlay::sequence<pargeo::point<5>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    return py_knnHelper<5>(array_vec, k, weighted);
  } else if (dim == 6) {
    parlay::sequence<pargeo::point<6>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    return py_knnHelper<6>(array_vec, k, weighted);
  } else if (dim == 7) {
    parlay::sequence<pargeo::point<7>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
    return py_knnHelper<7>(array_vec, k, weighted);
  } else
    throw std::runtime_error("Only dimensions 2-7 are supported");
}

PYBIND11_MODULE(pypargeo, m)
{
  m.doc() = "Pargeo: library for parallel computational geometry.";

  m.def("loadPoints", &py_loadPoints, "Load points based on file path.");

  m.def("KnnGraph", &py_knnGraph, "Input: 2d-numpy array containing n points in R^2; parameter k. Output: array E of size n*k, where E[i] and E[i+1] forall i are point indices represents a directed edge of the directed-knn graph.", py::arg("array"), py::arg("k"), py::arg("weighted")=false);

  m.def("DelaunayGraph", &py_delaunayGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of edges, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the delaunay graph.", py::arg("array"), py::arg("weighted")=false);

  m.def("GabrielGraph", &py_gabrielGraph, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the Gabriel graph.", py::arg("array"), py::arg("weighted")=false);

  m.def("BetaSkeleton", &py_skeleton, "Input: 2d-numpy array containing n points in R^2. Output: array E of size n*2, where E[i] and E[i+1] forall i are point indices that represents an undirected edge of the beta skeleton graph.", py::arg("array"), py::arg("beta"), py::arg("weighted")=false);
}
