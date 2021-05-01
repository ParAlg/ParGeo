project(kdtKnn)

add_executable(${PROJECT_NAME} knnTime.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)
