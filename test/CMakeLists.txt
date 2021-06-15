cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

# check whether googletest is locally installed, if not download and fetch
message(STATUS "--------------- GoogleTest -------------")
find_package(GTest CONFIG)
if(NOT GTest_FOUND)
    # new way of including googletest
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/googletest-download )
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/googletest-download )
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
            ${CMAKE_BINARY_DIR}/googletest-build
            EXCLUDE_FROM_ALL)
    set(GTest_LIBRARIES "gtest")
else()
    message(STATUS "using locally installed GoogleTest")
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(GTest_LIBRARIES GTest::gtest)
endif()

include(CTest)

add_executable(kdTree_test kdTree/kdTree_test.cpp)
target_link_libraries(kdTree_test PRIVATE
  pargeolib ${GTest_LIBRARIES})
add_test(NAME kdTree_test COMMAND kdTree_test)

add_executable(closestPair_test closestPair/closestPair_test.cpp)
target_link_libraries(closestPair_test PRIVATE
  pargeolib ${GTest_LIBRARIES})
add_test(NAME cloestPair_test COMMAND closestPair_test)

add_executable(wspd_test wspd/wspd_test.cpp)
target_link_libraries(wspd_test PRIVATE
  pargeolib ${GTest_LIBRARIES})
add_test(NAME wspd_test COMMAND wspd_test)

add_executable(hull_test hull/hull_test.cpp)
target_link_libraries(hull_test PRIVATE
  pargeolib convexHull3d ${GTest_LIBRARIES})
add_test(NAME hull_test COMMAND hull_test)

add_executable(hullTopology_test hull/hullTopology_test.cpp)
target_link_libraries(hullTopology_test PRIVATE
  pargeolib convexHull3d ${GTest_LIBRARIES})
add_test(NAME hullTopology_test COMMAND hullTest_test)

add_executable(dataset_test dataset/dataset_test.cpp)
target_link_libraries(dataset_test PRIVATE
  pargeolib dataset ${GTest_LIBRARIES})
add_test(NAME dataset_test COMMAND dataset_test)

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
file(COPY datasets DESTINATION ${CMAKE_BINARY_DIR}/test)
