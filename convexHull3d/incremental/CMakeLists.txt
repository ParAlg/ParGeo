set(CMAKE_CXX_FLAGS "-DRESERVE -DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")
add_executable(incremental hull.cpp hullTime.cpp)
