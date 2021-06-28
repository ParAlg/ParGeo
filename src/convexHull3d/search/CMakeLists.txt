add_library(searchHull3dLib hull.cpp)

target_link_libraries(searchHull3dLib pargeoLib bruteforceHull3dLib)

#target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
