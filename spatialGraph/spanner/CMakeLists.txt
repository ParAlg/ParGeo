project(euclideanSpanner)

add_executable(${PROJECT_NAME} euclideanSpanner.cpp graphTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
