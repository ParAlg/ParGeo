set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

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
# message(STATUS "benchmarks: Enabled")

add_executable(spatialGraph_bench spatialGraph_bench.cpp)
target_link_libraries(spatialGraph_bench PRIVATE
  pargeoLib
  spatialGraphLib
  emstLib
  benchmark_main)

add_executable(wspd_bench wspd_bench.cpp)
target_link_libraries(wspd_bench PRIVATE
  pargeoLib
  benchmark_main)

add_executable(zorder_bench zorder_bench.cpp)
target_link_libraries(zorder_bench PRIVATE
  pargeoLib
  benchmark_main)

add_executable(spatialSearch_bench spatialSearch_bench.cpp)
target_link_libraries(spatialSearch_bench PRIVATE
  pargeoLib
  benchmark_main)

add_executable(hull_bench hull_bench.cpp)
target_link_libraries(hull_bench PRIVATE
  divideConquerHull3dLib
  bruteforceHull3dLib
  samplingHull3dLib
  pseudoHull3dLib
  serialHull3dLib
  parallelHull3dLib
  pargeoLib
  benchmark_main)

add_executable(hull_sampling hull_sampling.cpp)
target_link_libraries(hull_sampling PRIVATE
  samplingHull3dLib
  pargeoLib
  benchmark_main)
