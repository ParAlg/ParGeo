project(uniformData)

add_executable(${PROJECT_NAME} generator.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
