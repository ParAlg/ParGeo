#include <string>
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

#include "dataset/seedSpreader.h"

template<int dim>
void ssGenerator(size_t n, size_t shape, char* fName) {
  if (shape == 0) {
    auto P = pargeo::seedSpreader::simdenGenerator<dim>(n);
    pargeo::pointIO::writePointsToFile(P, fName);
  } else {
    auto P = pargeo::seedSpreader::vardenGenerator<dim>(n);
    pargeo::pointIO::writePointsToFile(P, fName);
  }
}

int main(int argc, char* argv[]) {

  std::string text = "[-s {0: simden, 1: varden}] [-d {2--9}] n <outFile>";
  pargeo::commandLine P(argc, argv, text);

  std::pair<size_t, char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fName = in.second;

  int dim = P.getOptionIntValue("-d", 2);
  size_t shape = P.getOptionIntValue("-s", 0);

  if (dim == 2) ssGenerator<2>(n, shape, fName);
  else if (dim == 3) ssGenerator<3>(n, shape, fName);
  else if (dim == 4) ssGenerator<4>(n, shape, fName);
  else if (dim == 5) ssGenerator<5>(n, shape, fName);
  else if (dim == 6) ssGenerator<6>(n, shape, fName);
  else if (dim == 7) ssGenerator<7>(n, shape, fName);
  else if (dim == 8) ssGenerator<8>(n, shape, fName);
  else if (dim == 9) ssGenerator<9>(n, shape, fName);
  else throw std::runtime_error("dimensionality not yet supported");

  return 0;
}
