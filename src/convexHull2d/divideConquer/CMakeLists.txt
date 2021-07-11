add_library(divideHull2dLib
  hull.cpp
  )

target_link_libraries(divideHull2dLib
  pargeoLib
  bruteforceHull2dLib
  quickHull2dLib
  )
