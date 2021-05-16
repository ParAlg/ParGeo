project(wspd)

add_executable(${PROJECT_NAME} wspdTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)
