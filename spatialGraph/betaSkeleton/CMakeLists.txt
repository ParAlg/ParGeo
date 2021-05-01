project(betaSkeleton)

add_executable(${PROJECT_NAME} betaSkeleton.cpp delaunayGraph.cpp graphTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib delaunay2d)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
