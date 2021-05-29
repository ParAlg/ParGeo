cmake_minimum_required(VERSION 3.10)

project(clustering)

add_library(${PROJECT_NAME}
  hdbscan/hdbscan.cpp
  hdbscan/hdbscanTime.cpp
  dendrogram/dendrogram.cpp
  )

add_subdirectory(hdbscan)
add_subdirectory(dendrogram)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
