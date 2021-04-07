add_subdirectory(pybind11)
pybind11_add_module(pygraph pygraph.cpp knnGraph.cpp delaunayGraph.cpp gabrielGraph.cpp betaSkeleton.cpp incremental/delaunay.cpp)
