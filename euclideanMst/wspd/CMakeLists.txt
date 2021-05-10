project(emstWspd)

add_executable(${PROJECT_NAME} emstTime.cpp emst.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
