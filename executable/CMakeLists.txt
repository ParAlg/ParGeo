add_executable(closestPair closestPair.cpp)
target_link_libraries(closestPair pargeoLib)

add_executable(delaunay2d delaunay2d.cpp)
target_link_libraries(delaunay2d pargeoLib delaunay2dLib)

add_executable(graphGenerator graphGenerator.cpp)
target_link_libraries(graphGenerator
  pargeoLib
  delaunay2dLib
  emstLib
  spatialGraphLib)

add_executable(knn knn.cpp)
target_link_libraries(knn pargeoLib)

add_executable(dynKdTree dynKdTree.cpp)
target_link_libraries(dynKdTree pargeoLib)

add_executable(zorder zorder.cpp)
target_link_libraries(zorder pargeoLib)

add_executable(uniformData uniformData.cpp)
target_link_libraries(uniformData pargeoLib)

add_executable(hull hull.cpp)
target_link_libraries(hull pargeoLib
  divideHull2dLib
  pseudoHull3dLib
  )

add_executable(seb seb.cpp)
target_link_libraries(seb pargeoLib
  samplingSebLib
  welzlSebLib
  scanSebLib
  )
