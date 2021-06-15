cmake_minimum_required(VERSION 3.10)

project(pargeo_benchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# # User option to build benchmarks
# option(PARLAY_BENCHMARK    "Build microbenchmarks"  OFF)

# Benchmark should not run its own unit tests
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)

# Download Benchmark library
include(FetchContent)
FetchContent_Declare(benchmark
  GIT_REPOSITORY  https://github.com/google/benchmark
  GIT_TAG         master
)
FetchContent_GetProperties(benchmark)
if(NOT benchmark_POPULATED)
  message(STATUS "benchmarks: Configuring Google Benchmark")
  FetchContent_Populate(benchmark)
  set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE BOOL "")
  add_subdirectory(${benchmark_SOURCE_DIR}
                   ${benchmark_BINARY_DIR}
                   EXCLUDE_FROM_ALL)
endif()

# Include benchmark targets
message(STATUS "benchmarks: Enabled")

add_executable(spatialGraph_bench spatialGraph_bench.cpp)
target_link_libraries(spatialGraph_bench PRIVATE
  pargeolib
  dataset
  spatialGraph
  benchmark_main)

add_executable(hull_bench hull_bench.cpp)
target_link_libraries(hull_bench PRIVATE
  pargeolib
  dataset
  convexHull3d
  benchmark_main)
