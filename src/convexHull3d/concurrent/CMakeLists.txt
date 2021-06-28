add_library(concurrentHull3dLib
  hull.cpp
  )

target_link_libraries(concurrentHull3dLib
  pargeoLib
  serialHull3dLib
  incrementalHull3dLib
  )
