project(uniformData)

add_executable(${PROJECT_NAME} uniform.cpp generator.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
