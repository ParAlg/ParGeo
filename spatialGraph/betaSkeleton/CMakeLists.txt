project(betaSkeleton)

add_executable(${PROJECT_NAME} betaSkeleton.cpp graphTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib delaunay2d spatialGraph)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
