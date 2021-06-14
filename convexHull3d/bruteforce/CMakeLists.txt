project(bruteforceHull3d)

#------------------- Compile as executable

# for best serial performance uncomment below
# set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib serialHull3dLib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)

#------------------- Compile as library

add_library(${PROJECT_NAME}Lib
  hull.cpp
  )

target_link_libraries(${PROJECT_NAME}Lib pargeolib)

target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
