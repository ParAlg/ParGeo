project(octTreeHull3d)

set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -DGRID -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib serialHull3dLib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
