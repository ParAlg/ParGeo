#define VERBOSE

#ifdef VERBOSE
#include "common/get_time.h"
#endif

#include "kdtKnn.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<size_t> py_kdtknn(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k)
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
    parlay::sequence<size_t> result_vec = kdtKnn::kdtKnn<2, pargeo::point<2>>(array_vec, k);
#ifdef VERBOSE
    std::cout << "::knn-time = " << t.get_next() << std::endl;
#endif
    auto result = py::array_t<size_t>(k * n);
    auto result_buffer = result.request();
    size_t *result_ptr = (size_t *) result_buffer.ptr;
    std::memcpy(result_ptr,result_vec.data(),result_vec.size()*sizeof(size_t));
#ifdef VERBOSE
    std::cout << "::copy-out-time = " << t.stop() << std::endl;
#endif
    return result;
  } else if (dim == 3) {
#ifdef VERBOSE
    t.start();
#endif
    parlay::sequence<pargeo::point<3>> array_vec(n);
    std::memcpy(array_vec.data(), array.data(), array.size() * sizeof(double));
#ifdef VERBOSE
    std::cout << "::copy-in-time = " << t.get_next() << std::endl;
#endif
    parlay::sequence<size_t> result_vec = kdtKnn::kdtKnn<3, pargeo::point<3>>(array_vec, k);
#ifdef VERBOSE
    std::cout << "::knn-time = " << t.get_next() << std::endl;
#endif
    auto result = py::array_t<size_t>(k * n);
    auto result_buffer = result.request();
    size_t *result_ptr = (size_t *) result_buffer.ptr;
    std::memcpy(result_ptr,result_vec.data(),result_vec.size()*sizeof(size_t));
#ifdef VERBOSE
    std::cout << "::copy-out-time = " << t.stop() << std::endl;
#endif
    return result;
  } else
    throw std::runtime_error("Right now the code only supports 2d or 3d points.");

}

PYBIND11_MODULE(pyknn, m)
{
  m.doc() = "Kd-tree based Euclidean KNN on point data sets.";

  m.def("kdtknn", &py_kdtknn, "Input: 2d-numpy array containing points in R^2 or R^3; an integer k indicating the number of neighbors desired (k>=2 as the search includes self). Output: a 1d-numpy array of size k * n, containing indices of the KNNs of the input points concatenated.");
}
