cmake_minimum_required(VERSION 2.8.12)
project(example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(USE_PARALLEL "Parallel scheduler" ON)
option(USE_WARNINGS "Warning flags" OFF)

if (USE_WARNINGS)
  set(EXTRA_FLAGS "-Wall -Wextra")
endif()

if(USE_PARALLEL)
  set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")
endif()

find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_subdirectory(pybind11)
pybind11_add_module(example example.cpp)
