project(concurrentHull3d)

set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

#------------------- Compile as executable

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)

#------------------- Compile as library

add_library(${PROJECT_NAME}Lib
  hull.cpp
  )

target_link_libraries(${PROJECT_NAME}Lib pargeolib)

target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
