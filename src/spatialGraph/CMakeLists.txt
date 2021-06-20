add_library(spatialGraphLib
  delaunayGraph/delaunayGraph.cpp
  knnGraph/knnGraph.cpp
  gabrielGraph/gabrielGraph.cpp
  betaSkeleton/betaSkeleton.cpp
  spanner/spanner.cpp
  )

target_link_libraries(spatialGraphLib pargeoLib delaunay2dLib)
