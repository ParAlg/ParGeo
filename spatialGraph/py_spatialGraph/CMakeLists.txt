add_subdirectory(pybind11)
pybind11_add_module(geograph geograph.cpp knnGraph.cpp delaunayGraph.cpp gabrielGraph.cpp betaSkeleton.cpp incremental/delaunay.cpp)
