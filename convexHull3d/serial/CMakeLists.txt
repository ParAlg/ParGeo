project(serialHull3d)

#------------------- Compile as executable

# for best serial performance uncomment
# set(CMAKE_CXX_FLAGS "-DPARLAY_SEQUENTIAL -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

add_executable(${PROJECT_NAME} hullTime.cpp hull.cpp)

target_link_libraries(${PROJECT_NAME} pargeolib)

target_include_directories(${PROJECT_NAME} PRIVATE ../include)

#------------------- Compile as library

# set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

add_library(${PROJECT_NAME}Lib
  hull.cpp
  hullInternal.cpp
  hullInternal2.cpp
  hullInternal3.cpp
  )

target_link_libraries(${PROJECT_NAME}Lib pargeolib)

target_include_directories(${PROJECT_NAME}Lib PRIVATE ../include)
