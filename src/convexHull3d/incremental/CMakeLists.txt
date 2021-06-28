project(incrementalHull3d)

add_library(incrementalHull3dLib
  hull.cpp
  )

target_link_libraries(incrementalHull3dLib pargeoLib serialHull3dLib)

#target_include_directories(incrementalHull3dLib PRIVATE ../include)
