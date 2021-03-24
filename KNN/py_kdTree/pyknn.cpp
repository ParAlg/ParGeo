#include "kdtKnn.h"

template<int dim>
std::vector<size_t> knnCaller(std::vector<double> &P, size_t k) {
  if (k < 2 || k > P.size()) {
    std::cout << "Error, k = " << k << " k must range from 2 to #-data-points, abort" << std::endl;
    abort();
  }

  size_t n = P.size();

  double* mem = &P[0];
  auto memp = (pargeo::point<dim> *) mem;

  auto _P = parlay::sequence<pargeo::point<dim>>(n/dim);
  parlay::parallel_for(0, n/dim, [&](size_t i) {
			   _P[i] = memp[i];
			 });

  parlay::sequence<size_t> nnIdx = kdtKnn::kdtKnn<dim, pargeo::point<dim>>(_P, k);

  size_t *r = nnIdx.data();
  return std::vector(r, r+k*n/dim);
}

// ----------------
// Python interface
// ----------------
// reference: https://github.com/tdegeus/pybind11_examples/blob/master/03_numpy-1D_cpp-vector/example.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// wrap C++ function with NumPy array IO
py::array_t<size_t> py_kdtknn(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t k)
{
  // allocate std::vector (to pass to the C++ function)
  std::vector<double> array_vec(array.size());

  // copy py::array -> std::vector
  std::memcpy(array_vec.data(),array.data(),array.size()*sizeof(double));

  // call pure C++ function
  std::vector<size_t> result_vec = knnCaller<2>(array_vec, k); // todo dimensionality

  // allocate py::array (to pass the result of the C++ function to Python)
  auto result        = py::array_t<size_t>(array.size());
  auto result_buffer = result.request();
  size_t *result_ptr    = (size_t *) result_buffer.ptr;

  // copy std::vector -> py::array
  std::memcpy(result_ptr,result_vec.data(),result_vec.size()*sizeof(size_t));

  return result;
}

// wrap as Python module
PYBIND11_MODULE(pyknn, m)
{
  m.doc() = "pybind11 knn";

  m.def("kdtknn", &py_kdtknn, "Compute the knn of the input points");
}
