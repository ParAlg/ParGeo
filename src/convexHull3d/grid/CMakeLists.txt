add_library(gridHull3dLib
  hull.cpp
  )

target_link_libraries(gridHull3dLib
  pargeoLib
  serialHull3dLib
  concurrentHull3dLib
  )
