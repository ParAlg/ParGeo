cmake_minimum_required(VERSION 3.10)

project(euclideanMst)

add_library(${PROJECT_NAME}
  wspd/emst.cpp
  wspd/emstTime.cpp
  )

add_subdirectory(wspd)

add_subdirectory(filterKruskal)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
