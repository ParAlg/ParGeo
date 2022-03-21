cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)

file(GLOB SRCS *.cpp)
include(GoogleTest)
add_executable(test_binary_heap_layout ${SRCS})
target_link_libraries(test_binary_heap_layout PRIVATE
  pargeoLib
  ${GTest_LIBRARIES})

gtest_discover_tests(test_binary_heap_layout)
#add_test(NAME test_cache_oblivious)
