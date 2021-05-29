cmake_minimum_required(VERSION 3.10)

project(spatialSort)

add_library(${PROJECT_NAME}
  zorder/sortTime.cpp
  )

add_subdirectory(zorder)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
