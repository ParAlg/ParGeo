add_library(parallelHull3dLib
  hull.cpp
  )

target_link_libraries(parallelHull3dLib pargeoLib serialHull3dLib)
