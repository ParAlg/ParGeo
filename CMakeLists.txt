cmake_minimum_required(VERSION 3.10)

project(pargeo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(USE_PARALLEL "Parallel scheduler" ON)
option(USE_WARNINGS "Warning flags" OFF)

if (USE_WARNINGS)
  set(EXTRA_FLAGS "-Wall -Wextra")
endif()

if(USE_PARALLEL)
  set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG -fPIC ${EXTRA_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -mcx16 -O3 -DNDEBUG -fPIC ${EXTRA_FLAGS}")
endif()

message(STATUS "--------------- General configuration -------------")
message(STATUS "CMake Generator:                ${CMAKE_GENERATOR}")
message(STATUS "Compiler:                       ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build type:                     ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS:                ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG:          ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE:        ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS          ${CMAKE_CXX_LINKER_FLAGS}")
message(STATUS "CMAKE_INSTALL_PREFIX:           ${CMAKE_INSTALL_PREFIX}" )

# -------------------------------------------------
# Add benchmark directories

add_subdirectory(wspd)

add_subdirectory(knnSearch)

add_subdirectory(spatialSort)

add_subdirectory(clustering)

add_subdirectory(rangeSearch)

add_subdirectory(closestPair)

add_subdirectory(convexHull3d)

add_subdirectory(euclideanMst)

add_subdirectory(spatialGraph)

add_subdirectory(delaunayTriangulation)

add_subdirectory(pybindings)

add_subdirectory(dataset)

# -------------------------------------------------
# Testing

enable_testing()

add_subdirectory(test)

# -------------------------------------------------
# Benchmark

add_subdirectory(benchmark)

# -------------------------------------------------
# Link pargeolib as a library

project(pargeolib)

add_library(pargeolib INTERFACE)

target_include_directories(pargeolib INTERFACE include)
