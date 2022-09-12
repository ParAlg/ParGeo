set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Benchmark should not run its own unit tests
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)

# Download Benchmark library
include(FetchContent)
FetchContent_Declare(benchmark
  GIT_REPOSITORY  https://github.com/google/benchmark
  GIT_TAG         main
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
# message(STATUS "benchmarks: Enabled")

file(COPY datasets DESTINATION ${CMAKE_BINARY_DIR}/benchmark/)

add_executable(spatialGraph_bench spatialGraph_bench.cpp)
target_link_libraries(spatialGraph_bench PRIVATE
  pargeoLib
  spatialGraphLib
  emstLib
  benchmark_main)

add_executable(hull2d_bench hull2d_bench.cpp)
target_link_libraries(hull2d_bench PRIVATE
  randIncHull2dLib
  quickHull2dLib
  divideHull2dLib
  pargeoLib
  benchmark_main)

add_executable(hull3d_bench hull3d_bench.cpp)
target_link_libraries(hull3d_bench PRIVATE
  serialHull3dLib
  pseudoHull3dLib
  pargeoLib
  benchmark_main)

add_executable(seb_bench seb_bench.cpp)
target_link_libraries(seb_bench PRIVATE
  pargeoLib
  welzlSebLib
  scanSebLib
  samplingSebLib
  benchmark_main)

add_executable(all_bench all_bench.cpp)
target_link_libraries(all_bench PRIVATE
  emstLib
  divideHull2dLib
  pseudoHull3dLib
  samplingSebLib
  spatialGraphLib
  closestPairLib
  pargeoLib
  benchmark_main)
