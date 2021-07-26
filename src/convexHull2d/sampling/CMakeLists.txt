add_library(samplingHull2dLib
  hull.cpp
  )

target_link_libraries(samplingHull2dLib
  pargeoLib
  bruteforceHull2dLib
  quickHull2dLib
  )
