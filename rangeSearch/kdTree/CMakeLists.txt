project(kdtRangeSearch)

add_executable(${PROJECT_NAME} rangeTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)
