#pragma once

inline tuple<short,short,short> unpack(size_t g) {
  return make_tuple((short)g, (short)(g>>16), (short)(g>>32));
}

inline size_t pack(pargeo::fpoint<3> p, pargeo::fpoint<3>::floatT gsize) {
  size_t g = 0;
  for(int i=0; i<3; ++i) {
    unsigned short tmp = floor(p[i] / gsize);
    size_t tmp2 = 0;
    tmp2 += tmp;
    g |= tmp2 << (16*i);
  }
  return g;
}

struct gridAtt3d {
  size_t id;

  gridAtt3d() {}

  gridAtt3d(pargeo::fpoint<3> p, pargeo::fpoint<3>::floatT gsize) {
    id = pack(p, gsize);
  }
};

using gridpt3d = pargeo::_point<3, pointT::floatT, pointT::floatT, gridAtt3d>;
