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
# Add sub directories

add_subdirectory(src)

add_subdirectory(example)

add_subdirectory(benchmark)

enable_testing()

add_subdirectory(test)

# -------------------------------------------------
# Expose pargeoLib

project(pargeoLib)

add_library(pargeoLib INTERFACE)

target_include_directories(pargeoLib INTERFACE include)
