project(delaunay2d)

add_executable(delaunay incremental/delaunayTime.cpp incremental/delaunay.cpp)

add_library(${PROJECT_NAME} incremental/delaunay.cpp)

target_include_directories(delaunay
    PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_include_directories(${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(delaunay pargeolib)

target_link_libraries(${PROJECT_NAME} pargeolib)
