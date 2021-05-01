set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")
add_executable(serial hull.cpp hullTime.cpp)
