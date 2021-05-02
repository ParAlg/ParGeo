project(gabrielGraph)

add_executable(${PROJECT_NAME} gabrielGraph.cpp graphTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib delaunay2d)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
