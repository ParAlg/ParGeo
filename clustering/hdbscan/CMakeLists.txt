project(hdbscan)

add_executable(${PROJECT_NAME} hdbscanTime.cpp hdbscan.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib clustering)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
