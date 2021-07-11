add_executable(closestPair closestPair.cpp)
target_link_libraries(closestPair pargeoLib)

add_executable(delaunay2d delaunay2d.cpp)
target_link_libraries(delaunay2d pargeoLib delaunay2dLib)

add_executable(emst emst.cpp)
target_link_libraries(emst pargeoLib emstLib)

add_executable(betaSkeleton betaSkeleton.cpp)
target_link_libraries(betaSkeleton pargeoLib delaunay2dLib
  spatialGraphLib)

add_executable(delaunayGraph delaunayGraph.cpp)
target_link_libraries(delaunayGraph pargeoLib delaunay2dLib
  spatialGraphLib)

add_executable(gabrielGraph gabrielGraph.cpp)
target_link_libraries(gabrielGraph pargeoLib delaunay2dLib
  spatialGraphLib)

add_executable(knnGraph knnGraph.cpp)
target_link_libraries(knnGraph pargeoLib spatialGraphLib)

add_executable(spanner spanner.cpp)
target_link_libraries(spanner pargeoLib spatialGraphLib)

add_executable(knn knn.cpp)
target_link_libraries(knn pargeoLib)

add_executable(zorder zorder.cpp)
target_link_libraries(zorder pargeoLib)

add_executable(uniformData uniformData.cpp)
target_link_libraries(uniformData pargeoLib)

add_executable(hull hull.cpp)
target_link_libraries(hull pargeoLib
  divideConquerHull3dLib
  bruteforceHull3dLib
  samplingHull3dLib
  pseudoHull3dLib
  serialHull3dLib
  parallelHull3dLib
  )

add_executable(hull2d hull2d.cpp)
target_link_libraries(hull2d pargeoLib
  quickHull2dLib
  bruteforceHull2dLib
  divideHull2dLib
  )

add_executable(hullSampling hullSampling.cpp)
target_link_libraries(hullSampling pargeoLib
  samplingHull3dLib
  serialHull3dLib
  )
