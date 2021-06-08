project(bruteforceHull3d)

#------------------- Compile as executable

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib serialHull3dLib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)
