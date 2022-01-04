add_library(pseudoHull3dLib
  hull.cpp
  )

target_link_libraries(pseudoHull3dLib
  pargeoLib
  parallelHull3dLib
  )
