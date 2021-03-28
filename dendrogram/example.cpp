#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <limits> // numeric_limits
#include <algorithm> // nth_element
#include "common/pbbs/gettime.h"
#include "unionFind.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"


namespace py = pybind11;

// namespace dendrogram{
struct edge{
  pair<int, int> e;
  double w;
  int n=0;
  edge(int a, int b, double c):w(c){
    e = make_pair(a,b);
  }
  int first(){return e.first;}
  int second(){return e.second;}
  void set(int a, int b, int s){
     e = make_pair(a,b);
     n = s;
  }
  edge(){}
};

py::array_t<double> dendrogram(py::array_t<int, py::array::c_style | py::array::forcecast> edges, py::array_t<double, py::array::c_style | py::array::forcecast> weights)
{

  if (edges.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

  if(edges.shape()[1] != 2)
    throw std::runtime_error("Input should be edge lists");

  size_t n = (edges.size() / 2) + 1; //number of vertices

  unionFind uf = unionFind(n);

  int idx = n;
  parlay::sequence<int> idxmap(n);
  parlay::sequence<int> sizes(n);
  parlay::parallel_for(0, n,[&](int i) {idxmap[i] = i;sizes[i] = 1;});

  parlay::sequence<edge> edges_sorted(n-1);
  auto ptr = static_cast<int *>(edges.request().ptr);
  auto ptr2 = static_cast<double *>(weights.request().ptr);
  parlay::parallel_for(0, n-1,
    [&](int i) {
      if(*(ptr+2*i) <  *(ptr+2*i+1)){
        edges_sorted[i] = edge(*(ptr+2*i), *(ptr+2*i+1), *(ptr2+i));}
      else{
        edges_sorted[i] = edge(*(ptr+2*i+1), *(ptr+2*i), *(ptr2+i));}
  });

  parlay::internal::sample_sort_inplace(parlay::make_slice(edges_sorted), [&](edge i, edge j){return i.w < j.w;});

  for(int i=0; i<n-1; ++i){
    int u = uf.find(edges_sorted[i].first());
    int v = uf.find(edges_sorted[i].second());
    edges_sorted[i].set(idxmap[u], idxmap[v], sizes[u]+sizes[v]);
    uf.link(u, v);
    int newc = uf.find(u);
    idxmap[newc] = idx;
    sizes[newc] = sizes[u]+sizes[v];
    idx++;
  }
  uf.del();
  idxmap.clear();
  sizes.clear();

  auto result = py::array_t<double>(4 * (n-1));
  auto result_buffer = result.request();
  double *result_ptr = (double *) result_buffer.ptr;
  parlay::parallel_for(0, n-1,
    [&](int i) {
      *(result_ptr+4*i) = static_cast<double>(edges_sorted[i].first());
      *(result_ptr+4*i+1) = static_cast<double>(edges_sorted[i].second());
      *(result_ptr+4*i+2) = static_cast<double>(edges_sorted[i].w);
      *(result_ptr+4*i+3) = static_cast<double>(edges_sorted[i].n);
  });

  edges_sorted.clear();

#ifdef VERBOSE
  std::cout << "::copy-out-time = " << t.stop() << std::endl;
#endif
  return result;

}

int add(int i, int j) {
    return i + j;
}


class Pet
{
    public:
        Pet(const std::string &name, int hunger) : name(name), hunger(hunger) {}
        ~Pet() {}

        void go_for_a_walk() { hunger++; }
        const std::string &get_name() const { return name; }
        int get_hunger() const { return hunger; }

    private:
        std::string name;
        int hunger;
};

// }

PYBIND11_MODULE(example, m) {
    // optional module docstring
    m.doc() = "pybind11 example plugin";

    // define add function
    m.def("add", &add, "A function which adds two numbers");
    m.def("dendrogram", &dendrogram, "A function which produce a dendrogram");

    // bindings to Pet class
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &, int>())
        .def("go_for_a_walk", &Pet::go_for_a_walk)
        .def("get_hunger", &Pet::get_hunger)
        .def("get_name", &Pet::get_name);
}
