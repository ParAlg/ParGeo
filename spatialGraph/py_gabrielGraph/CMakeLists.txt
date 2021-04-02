add_subdirectory(pybind11)
pybind11_add_module(pygabriel pygabriel.cpp gabrielGraph.cpp incremental/delaunay.cpp)
