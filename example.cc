#include "qgb.hh"

#include <algorithm>
#include <fstream>

using namespace Geometry;

QGB::Boundary readCurve(std::istream &is) {
  QGB::Boundary b;
  for (size_t i = 0; i <= 2; ++i)
    is >> b[i];
  return b;
}

QGB readPatch(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n;
  f >> n;
  QGB result(n);

  for (size_t i = 0; i < n; ++i)
    result.setBoundary(i, readCurve(f));

  try {
    Point3D mp;
    f >> mp;
    result.setMidpoint(mp);
  } catch(std::ios_base::failure &) {
    result.resetMidpoint();
  }

  return result;
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.qgb> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  auto patch = readPatch(argv[1]);
  patch.eval(resolution).writeOBJ("test.obj");
}
