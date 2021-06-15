cmake_minimum_required(VERSION 3.10)

project(dataset)

add_library(${PROJECT_NAME}
  uniform/uniform.cpp
  )

add_subdirectory(uniform)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)
