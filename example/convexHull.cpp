#include "convexHull2d/divideConquer/hull.h"
#include "convexHull3d/pseudo/hull.h"

#include <iostream>
#include <algorithm>
#include "pargeo/parseCommandLine.h"
#include "parlay/parallel.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include "pargeo/pointIO.h"

using namespace pargeo;
using namespace pargeo::pointIO;

int main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-o <outFile>] <inFile>\n hull2d output format: point indices\n hull3d output format: #facets * 9 floats representing coordinates of facet vertices");
  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");

  int dim = readHeader(iFile);

  switch(dim) {
  case 2:
    {
      parlay::sequence<pargeo::fpoint<2>> P = readPointsFromFile<pargeo::fpoint<2>>(iFile);
      auto H = pargeo::hull2d::divideConquer::compute<pargeo::fpoint<2>>(make_slice(P));
      pargeo::IO::writeSeqToFile("", H, oFile);
      break;
    }
  case 3:
    {
      parlay::sequence<pargeo::fpoint<3>> P = readPointsFromFile<pargeo::fpoint<3>>(iFile);
      auto H = pargeo::hull3d::pseudo::compute<pargeo::fpoint<3>>(make_slice(P));

      if (oFile != NULL) {

	parlay::sequence<double> pOut(H.size() * 3 * 3);

	parlay::parallel_for(0, H.size(), [&](size_t i) {
					    for (int d = 0; d < 3; ++ d)
					      pOut[i * 9 + d] = H[i].a[d];

					    for (int d = 0; d < 3; ++ d)
					      pOut[i * 9 + 3 + d] = H[i].b[d];

					    for (int d = 0; d < 3; ++ d)
					      pOut[i * 9 + 6 + d] = H[i].c[d];
					  });

	pargeo::pointIO::writePointsToFile(pOut, oFile);

      }

      break;
    }
  default:
    throw std::runtime_error("Convex hull only supports 2d and 3d data sets.");
  }

}
