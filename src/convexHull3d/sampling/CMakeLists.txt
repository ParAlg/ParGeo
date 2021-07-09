add_library(samplingHull3dLib
  hull.cpp
  hullTest.cpp)

target_link_libraries(samplingHull3dLib
  serialHull3dLib
  parallelHull3dLib
  pargeoLib)
