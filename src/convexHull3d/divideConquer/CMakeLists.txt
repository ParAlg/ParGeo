add_library(divideConquerHull3dLib
  hull.cpp
  )

target_link_libraries(divideConquerHull3dLib
  pargeoLib
  samplingHull3dLib
  serialHull3dLib
  parallelHull3dLib
  )
