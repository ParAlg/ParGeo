#include <string>
#include "pargeo/point.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

#include "dataset/uniform.h"

template<int dim>
void uniformGenerator(size_t n, size_t shape, double thickness, char* fName) {
  if (thickness < 0) {
    auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, shape, sqrt(double(n)));
    pargeo::pointIO::writePointsToFile(P, fName);
  } else {
    auto P = pargeo::uniformOnPolyPoints<dim, pargeo::point<dim>>(n, shape, thickness, double(n));
    pargeo::pointIO::writePointsToFile(P, fName);
  }
}

int main(int argc, char* argv[]) {
  //using namespace pargeo;

  std::string text = "[-s] [-c] [-t {double}] [-d {2--9}] n <outFile>";
  text += "\n polygon type: -s sphere -c cube";
  text += "\n generate point on surface: -t thickness";
  text += "\n  o.w. default to generate in polygon";
  pargeo::commandLine P(argc, argv, text);

  std::pair<size_t, char*> in = P.sizeAndFileName();
  size_t n = in.first;
  char* fName = in.second;

  int dim = P.getOptionIntValue("-d", 2);
  bool sphere = P.getOption("-s");
  bool cube = P.getOption("-c");
  double thickness = P.getOptionDoubleValue("-t", -1);

  size_t shape = sphere ? 0 : 1;

  if (dim == 2) uniformGenerator<2>(n, shape, thickness, fName);
  else if (dim == 3) uniformGenerator<3>(n, shape, thickness, fName);
  else if (dim == 4) uniformGenerator<4>(n, shape, thickness, fName);
  else if (dim == 5) uniformGenerator<5>(n, shape, thickness, fName);
  else if (dim == 6) uniformGenerator<6>(n, shape, thickness, fName);
  else if (dim == 7) uniformGenerator<7>(n, shape, thickness, fName);
  else if (dim == 8) uniformGenerator<8>(n, shape, thickness, fName);
  else if (dim == 9) uniformGenerator<9>(n, shape, thickness, fName);
  else throw std::runtime_error("dimensionality not yet supported");

  return 0;
}
