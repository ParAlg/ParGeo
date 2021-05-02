cmake_minimum_required(VERSION 3.10)

project(spatialGraph)

add_library(${PROJECT_NAME}
  delaunayGraph/delaunayGraph.cpp
  knnGraph/knnGraph.cpp
  gabrielGraph/gabrielGraph.cpp
  betaSkeleton/betaSkeleton.cpp
  )

add_subdirectory(knnGraph)
add_subdirectory(delaunayGraph)
add_subdirectory(gabrielGraph)
add_subdirectory(betaSkeleton)

target_link_libraries(${PROJECT_NAME} pargeolib delaunay2d)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
