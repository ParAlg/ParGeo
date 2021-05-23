project(pypargeo)

add_subdirectory(pybind11)

pybind11_add_module(${PROJECT_NAME} pypargeo.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE pargeolib spatialGraph delaunay2d clustering)
