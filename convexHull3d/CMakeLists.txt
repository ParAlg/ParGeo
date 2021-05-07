cmake_minimum_required(VERSION 3.10)

project(convexHull3d)

add_library(${PROJECT_NAME}
  concurrent/hull.cpp
  concurrent/hullTime.cpp
  giftWrap/hullTime.cpp
  octTree/hull.cpp
  octTree/hullTime.cpp
  )

add_subdirectory(serial)
add_subdirectory(incremental)
add_subdirectory(concurrent)
add_subdirectory(octTree)
add_subdirectory(giftWrap)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
