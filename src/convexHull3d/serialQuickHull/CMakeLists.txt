# set(CMAKE_CXX_FLAGS "-DHOMEGROWN -pthread -mcx16 -O3 -DNDEBUG ${EXTRA_FLAGS}")

add_library(serialHull3dLib
  hull.cpp
  )

target_link_libraries(serialHull3dLib
  pargeoLib
  bruteforceHull3dLib
  )

target_include_directories(serialHull3dLib PRIVATE header)
