cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

file(GLOB SRCS *.cpp)
include(GoogleTest)
add_executable(kdTree_test ${SRCS})
target_link_libraries(kdTree_test PRIVATE
  pargeolib
  ${GTest_LIBRARIES})
