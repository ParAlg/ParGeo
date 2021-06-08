cmake_minimum_required(VERSION 3.10)

project(convexHull3d)

add_library(${PROJECT_NAME}
  serial/hull.cpp
  serial/hullInternal.cpp
  serial/hullInternal2.cpp
  serial/hullInternal3.cpp
  octTree/hull.cpp
  concurrent/hull.cpp
  bruteforce/hull.cpp
  pseudoHull/hull.cpp
  incremental/hull.cpp
  sampling/hull.cpp
  sampling/hullTest.cpp
  )

add_subdirectory(incremental)
add_subdirectory(bruteforce)
add_subdirectory(pseudoHull)
add_subdirectory(concurrent)
add_subdirectory(sampling)
add_subdirectory(giftWrap)
add_subdirectory(octTree)
add_subdirectory(serial)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
  PUBLIC ${PROJECT_SOURCE_DIR}/include)
