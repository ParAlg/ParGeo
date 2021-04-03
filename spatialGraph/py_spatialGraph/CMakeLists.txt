add_subdirectory(pybind11)
pybind11_add_module(pygraph pygraph.cpp delaunayGraph.cpp gabrielGraph.cpp incremental/delaunay.cpp)
