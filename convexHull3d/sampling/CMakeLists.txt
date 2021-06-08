project(samplingHull3d)

#------------------- Compile as executable

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib serialHull3dLib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)

#------------------- Compile as library

add_library(${PROJECT_NAME}Lib
  hull.cpp
  hullTest.cpp
  )

target_link_libraries(${PROJECT_NAME}Lib pargeolib serialHull3dLib)

target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
