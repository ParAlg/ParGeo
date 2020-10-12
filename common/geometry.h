#ifndef _BENCH_GEOM_INCLUDED
#define _BENCH_GEOM_INCLUDED

#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <limits>
#include "pbbs/parallel.h"

using namespace std;

// *************************************************************
//    any dimensional POINTS
// *************************************************************

template <int _dim> class point;

template <int _dim> class vect {
public:
  typedef double floatT;
  typedef vect vectT;
  typedef point<_dim> pointT;
  floatT x[_dim];
  static const int dim = _dim;
  vect() { for (int i=0; i<_dim; ++i) x[i] = 0; }
  vect(pointT p) { for (int i=0; i<_dim; ++i) x[i] = p.x[i]; }
  vect(floatT* v) { for (int i=0; i<_dim; ++i) x[i] = v[i]; }
  void print() {
    cout << std::setprecision(2);
    cout << ":(";
    for (int i=0; i<_dim-1; ++i) cout << x[i] << ",";
    cout << x[_dim-1] << "):";
  }
  vectT operator+(vectT op2) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i]=x[i] + op2.x[i];
    return vectT(xx);
  }
  vectT operator-(vectT op2) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i]=x[i] - op2.x[i];
    return vectT(xx);
  }
  pointT operator+(pointT op2) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i]=x[i] + op2.x[i];
    return point<dim>(xx);
  }
  vectT operator*(floatT s) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i] = x[i]*s;
    return vectT(xx);
  }
  vectT operator/(floatT s) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i] = x[i]/s;
    return vectT(xx);
  }
  floatT& operator[] (int i) {return x[i];}
  float dot(vectT v) {
    floatT xx=0;
    for (int i=0; i<_dim; ++i) xx += x[i]*v[i];
    return xx;
  }
  vectT cross(vectT v) {
    cout << "lack of standard implementation of cross product fo R^d, see point and vector in 2d and 3d, abort" << endl; abort();
  }
  floatT maxDim() {
    floatT xx = x[0];
    for (int i=1; i<_dim; ++i) xx = max(xx,x[i]);
    return xx;
  }
  floatT length() {
    floatT xx=0;
    for (int i=0; i<_dim; ++i) xx += x[i]*x[i];
    return sqrt(xx);
  }
};

template <int _dim> class point {
public:
  typedef double floatT;
  typedef vect<_dim> vectT;
  typedef point pointT;
  floatT x[_dim];
  static const int dim = _dim;
  static constexpr double empty = numeric_limits<double>::max();
  int dimension() {return _dim;}
  void setEmpty() {x[0]=empty;}
  bool isEmpty() {return x[0]==empty;}
  point() { for (int i=0; i<_dim; ++i) x[i]=empty; }
  point(vectT v) { for (int i=0; i<_dim; ++i) x[i]=v[i]; }
  point(floatT* p) { for (int i=0; i<_dim; ++i) x[i]=p[i]; }
  point(pointT* p) { for (int i=0; i<_dim; ++i) x[i]=p->x[i]; }
  void print() {//Deprecate
    cout << std::setprecision(2);
    cout << ":(";
    for (int i=0; i<_dim-1; ++i) cout << x[i] << ",";
    cout << x[_dim-1] << "):";
  }
  pointT operator+(vectT op2) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i] = x[i]+op2.x[i];
    return pointT(xx);
  }
  pointT operator-(pointT op2) {
    floatT xx[_dim];
    for (int i=0; i<_dim; ++i) xx[i] = x[i]-op2.x[i];
    return pointT(xx);
  }
  floatT& operator[](int i) {return x[i];}
  friend bool operator==(pointT a, pointT b) {
    for (intT ii=0; ii<dim; ++ii) {
      if (abs(a[ii]-b[ii]) > 0.0) return false;}
    return true;
  }
  friend bool operator!=(pointT a, pointT b) {return !(a==b);}
  floatT* coordinate() {return x;}
  floatT coordinate(int i) {return x[i];}
  void updateX(int i, floatT val) {x[i]=val;}//Deprecate
  void updateCoordinate(int i, floatT val) {x[i]=val;}
  pointT average(pointT p2) {
    auto pp = pointT();
    for (int i=0; i<_dim; ++i) pp.x[i] = (p2[i] + x[i])/2;
    return pp;
  }
  void minCoords(pointT b) {
    for (int i=0; i<_dim; ++i) x[i] = min(x[i], b.x[i]);}
  void minCoords(floatT* b) {
    for (int i=0; i<_dim; ++i) x[i] = min(x[i], b[i]);}
  void maxCoords(pointT b) {
    for (int i=0; i<_dim; ++i) x[i] = max(x[i], b.x[i]);}
  void maxCoords(floatT* b) {
    for (int i=0; i<_dim; ++i) x[i] = max(x[i], b[i]);}
  intT quadrant(pointT center) {
    intT index = 0;
    intT offset = 1;
    for (int i=0; i<_dim; ++i) {
      if (x[i] > center.x[i]) index += offset;
      offset *= 2;
    }
    return index;
  }
  bool outOfBox(pointT center, floatT hsize) {
    for (int i=0; i<_dim; ++i) {
      if (x[i]-hsize > center.x[i] || x[i]+hsize < center.x[i])
        return true;
    }
    return false;
  }
  floatT pointDist(pointT p) {//deprecate
    floatT xx=0;
    for (int i=0; i<_dim; ++i) xx += (x[i]-p.x[i])*(x[i]-p.x[i]);
    return sqrt(xx);
  }
  floatT dist(pointT p) {
    floatT xx=0;
    for (int i=0; i<_dim; ++i) xx += (x[i]-p.x[i])*(x[i]-p.x[i]);
    return sqrt(xx);
  }
  floatT dot(pointT p2) {
    floatT r = 0;
    for(int i=0; i<dim; ++i) r += x[i]*p2[i];
    return r;}
  pointT mult(floatT c) {
    pointT r;
    for(int i=0; i<dim; ++i) r[i] = x[i]*c;
    return r;}
  pointT normalize() {
    auto r = pointT(x);
    floatT s = 0;
    for (int i=0; i<dim; ++i) s += x[i]*x[i];
    s = sqrt(s);
    for (int i=0; i<dim; ++i) r[i] /= s;
    return r;
  }
};

template <int dim>
static std::ostream& operator<<(std::ostream& os, const vect<dim> v) {
  for (int i=0; i<dim; ++i)
    os << v.x[i] << " ";
  return os;
}

template <int dim>
static std::ostream& operator<<(std::ostream& os, const point<dim> v) {
  for (int i=0; i<dim; ++i)
    os << v.x[i] << " ";
  return os;
}

// *************************************************************
//    specialized POINTS AND VECTORS (2d, 3d)
// *************************************************************

class point2d;

class vect2d {
public: 
  typedef double floatT;
  typedef vect<2> vectT;
  vectT v;
  inline floatT x() {return v.x[0];}
  inline floatT y() {return v.x[1];}
  vect2d(floatT xx,floatT yy) {v.x[0]=xx; v.x[1]=yy;}
  vect2d() {
    v = vect<2>();
    // v.x[0] = 0;
    // v.x[1] = 0;
  }
  vect2d(point2d p);
  vect2d(vectT vv) : v(vv) {};
  vect2d(floatT* p) {v.x[0]=p[0]; v.x[1]=p[1];};
  vect2d operator+(vect2d op2) {return vect2d(v+op2.v);}
  vect2d operator-(vect2d op2) {return vect2d(v-op2.v);}
  point2d operator+(point2d op2);
  vect2d operator*(floatT s) {return vect2d(v*s);}
  vect2d operator/(floatT s) {return vect2d(v/s);}
  floatT& operator[] (int i) {return v[i];};
  floatT dot(vect2d vv) {return v.dot(vv.v);}
  floatT cross(vect2d vv) {return x()*vv.y() - y()*vv.x();;}
  floatT maxDim() {return v.maxDim();}
  void print() {v.print();}
  floatT Length(void) { return v.length();}
  static const int dim = 2;
};

class point2d {
public: 
  typedef double floatT;
  typedef point<2> pointT;
  typedef vect2d vectT;
  //floatT x; floatT y;
  pointT p;
  inline floatT x() {return p.x[0];};
  inline floatT y() {return p.x[1];};
  int dimension() {return 2;}//Deprecate
  point2d(floatT xx,floatT yy) {
    p.x[0] = xx;
    p.x[1] = yy;}
  point2d() {
    p = pointT();
    // p.x[0] = 0;
    // p.x[1] = 0;
  }
  point2d(pointT pp): p(pp) {};
  point2d(vect2d v) {p.x[0] = v.x(); p.x[1] = v.y();};
  point2d(floatT* pp) {p.x[0] = pp[0]; p.x[1] = pp[1];};
  void print() {cout << ":(" << x() << "," << y() << "):";}
  //todo s
  point2d operator+(vect2d op2) {return point2d(x() + op2.x(), y() + op2.y());}//todo
  vect2d operator-(point2d op2) {return vect2d(x() - op2.x(), y() - op2.y());}
  floatT& operator[] (int i) {return p[i];};
  point2d minCoords(point2d b) {
    return point2d(min(x(),b.x()),min(y(),b.y()));}
  point2d maxCoords(point2d b) {
    return point2d(max(x(),b.x()),max(y(),b.y()));}
  int quadrant(point2d center) {
    return p.quadrant(center.p);//check
    // int index = 0;
    // if (x() > center.x()) index += 1;
    // if (y() > center.y()) index += 2;
    // return index;
  }
  // returns a point2d offset by offset in one of 4 directions 
  // depending on dir (an integer from [0..3])
  point2d offsetPoint(int dir, floatT offset) {
    floatT xx = x() + ((dir & 1) ? offset : -offset);
    floatT yy = y() + ((dir & 2) ? offset : -offset);
    return point2d(xx,yy);
  }
  bool outOfBox(point2d pt, floatT hsize) {
    return p.outOfBox(pt.p, hsize);
    // return ((x() - hsize > pt.x()) || (x() + hsize < pt.x()) ||
    //         (y() - hsize > pt.y()) || (y() + hsize < pt.y()));
  }
  floatT pointDist(point2d p2) {
    return p.pointDist(p2.p);
    //return sqrt((p.x-x)*(p.x-x)+(p.y-y)*(p.y-y));
  }
  static const int dim = 2;
};

inline vect2d::vect2d(point2d p) {v=vectT(p.p);}

inline point2d vect2d::operator+(point2d op2) {
  auto tmp = v+op2.p;
  return point2d(tmp.x);}

static std::ostream& operator<<(std::ostream& os, vect2d v) {
  return os << v.x() << " " << v.y();}

static std::ostream& operator<<(std::ostream& os, point2d v) {
  return os << v.x() << " " << v.y(); }

class point3d;

class vect3d {
public: 
  typedef double floatT;
  typedef vect<3> vectT;
  typedef point<3> pointT;
  //floatT x; floatT y; floatT z;
  vectT v;
  inline floatT x() {return v[0];};
  inline floatT y() {return v[1];};
  inline floatT z() {return v[2];};
  vect3d(floatT xx,floatT yy, floatT zz) {
    v.x[0] = xx; v.x[1] = yy; v.x[2] = zz;}
  vect3d() {v = vectT();}
  vect3d(point3d p);
  vect3d(vectT vv): v(vv) {};
  vect3d(pointT pp) {v = vectT(pp.x);};
  vect3d(floatT* p) {
    v.x[0] = p[0]; v.x[1] = p[1]; v.x[2] = p[2];}
  vect3d operator+(vect3d op2) {
    return vect3d(v+op2.v);}
  vect3d operator-(vect3d op2) {
    return vect3d(v-op2.v);}
  point3d operator+(point3d op2);
  vect3d operator*(floatT s) {return vect3d(v*s);}
  vect3d operator/(floatT s) {return vect3d(v/s);}
  floatT& operator[] (int i) {return v[i];}
  floatT dot(vect3d vv) {return v.dot(vv.v);}
  vect3d cross(vect3d vv) {
    return vect3d(y()*vv.z() - z()*vv.y(), z()*vv.x() - x()*vv.z(), x()*vv.y() - y()*vv.x());
  }
  floatT maxDim() {return v.maxDim();}
  void print() {cout << std::setprecision(10) << ":(" << x() << "," << y() << "," << z() << "):";}
  floatT Length(void) {return v.length();}
  static const int dim = 3;
};

class point3d {
public: 
  typedef double floatT;
  typedef point<3> pointT;
  typedef vect3d vectT;
  //floatT x; floatT y; floatT z;
  pointT p;
  inline floatT x() {return p[0];}
  inline floatT y() {return p[1];}
  inline floatT z() {return p[2];}
  int dimension() {return 3;}
  point3d(floatT xx,floatT yy, floatT zz) {
    p.x[0] = xx; p.x[1] = yy; p.x[2] = zz;}
  point3d() {p = pointT();}
  point3d(vect3d vv) {
    p = pointT(vv.v);
  };
  point3d(floatT* pp) {
    p.x[0] = pp[0]; p.x[1] = pp[1]; p.x[2] = pp[2];
  };
  void print() {cout << ":(" << x() << "," << y() << "," << z() << "):";}
  vect3d operator-(point3d op2) {
    return vect3d(p-op2.p);}
  point3d operator+(vect3d op2) {
    return point3d(p+op2.v);}
  point3d minCoords(point3d b) {
    return point3d(min(x(),b.x()),min(y(),b.y()),min(z(),b.z()));}
  point3d maxCoords(point3d b) {
    return point3d(max(x(),b.x()),max(y(),b.y()),min(z(),b.z()));}
  floatT& operator[] (int i) {return p[i];}
  int quadrant(point3d center) {
    // int index = 0;
    // if (x > center.x) index += 1;
    // if (y > center.y) index += 2;
    // if (z > center.z) index += 4;
    // return index;
    return p.quadrant(center.p);
  }
  // returns a point3d offset by offset in one of 8 directions 
  // depending on dir (an integer from [0..7])
  point3d offsetPoint(int dir, floatT offset) {
    floatT xx = x() + ((dir & 1) ? offset : -offset);
    floatT yy = y() + ((dir & 2) ? offset : -offset);
    floatT zz = z() + ((dir & 4) ? offset : -offset);
    return point3d(xx,yy,zz);
  }
  // checks if pt is outside of a box centered at this point with
  // radius hsize
  bool outOfBox(point3d pt, floatT hsize) { 
    // return ((x - hsize > pt.x) || (x + hsize < pt.x) ||
    //         (y - hsize > pt.y) || (y + hsize < pt.y) ||
    //         (z - hsize > pt.z) || (z + hsize < pt.z));
    return p.outOfBox(pt.p, hsize);
  }
  static const int dim = 3;
};

inline vect3d::vect3d(point3d p) {v=vectT(p.p);}

inline point3d vect3d::operator+(point3d op3) {
  auto tmp = v+op3.p;
  return point3d(tmp.x);}

//inline vect3d::vect3d(point3d p) { x = p.x; y = p.y; z = p.z;}

static std::ostream& operator<<(std::ostream& os, vect3d v) {
  return os << v.x() << " " << v.y() << " " << v.z(); }

static std::ostream& operator<<(std::ostream& os, point3d v) {
  return os << v.x() << " " << v.y() << " " << v.z();
}

// *************************************************************
//    GEOMETRY
// *************************************************************

// Returns twice the area of the oriented triangle (a, b, c)
inline double triArea(point2d a, point2d b, point2d c) {
  return (b-a).cross(c-a);
}

inline double triAreaNormalized(point2d a, point2d b, point2d c) {
  return triArea(a,b,c)/((b-a).Length()*(c-a).Length());
}

// Returns TRUE if the points a, b, c are in a counterclockise order
inline int counterClockwise(point2d a, point2d b, point2d c) {
  return (b-a).cross(c-a) > 0.0;
}

//takes row-major 3by3 matrix m
inline floatT determinant3by3(floatT* m) {
  return m[0*3+0] * (m[1*3+1] * m[2*3+2] - m[2*3+1] * m[1*3+2]) -
    m[0*3+1] * (m[1*3+0] * m[2*3+2] - m[1*3+2] * m[2*3+0]) +
    m[0*3+2] * (m[1*3+0] * m[2*3+1] - m[1*3+1] * m[2*3+0]);
}

//takes row-major 3by3 matrix m, output replaces m
inline void inverse3by3(floatT* m) {
  floatT invDet = 1/determinant3by3(m);
  floatT mInv[9];
  mInv[0*3+0] = (m[1*3+1] * m[2*3+2] - m[2*3+1] * m[1*3+2]) * invDet;
  mInv[0*3+1] = (m[0*3+2] * m[2*3+1] - m[0*3+1] * m[2*3+2]) * invDet;
  mInv[0*3+2] = (m[0*3+1] * m[1*3+2] - m[0*3+2] * m[1*3+1]) * invDet;
  mInv[1*3+0] = (m[1*3+2] * m[2*3+0] - m[1*3+0] * m[2*3+2]) * invDet;
  mInv[1*3+1] = (m[0*3+0] * m[2*3+2] - m[0*3+2] * m[2*3+0]) * invDet;
  mInv[1*3+2] = (m[1*3+0] * m[0*3+2] - m[0*3+0] * m[1*3+2]) * invDet;
  mInv[2*3+0] = (m[1*3+0] * m[2*3+1] - m[2*3+0] * m[1*3+1]) * invDet;
  mInv[2*3+1] = (m[2*3+0] * m[0*3+1] - m[0*3+0] * m[2*3+1]) * invDet;
  mInv[2*3+2] = (m[0*3+0] * m[1*3+1] - m[1*3+0] * m[0*3+1]) * invDet;
  for(int i=0; i<9; ++i) m[i] = mInv[i];
}

//takes row-major 4by4 matrix m
inline floatT determinant4by4(floatT* m) {
  return
    m[0*4+3] * m[1*4+2] * m[2*4+1] * m[3*4+0] - m[0*4+2] * m[1*4+3] * m[2*4+1] * m[3*4+0] -
    m[0*4+3] * m[1*4+1] * m[2*4+2] * m[3*4+0] + m[0*4+1] * m[1*4+3] * m[2*4+2] * m[3*4+0] +
    m[0*4+2] * m[1*4+1] * m[2*4+3] * m[3*4+0] - m[0*4+1] * m[1*4+2] * m[2*4+3] * m[3*4+0] -
    m[0*4+3] * m[1*4+2] * m[2*4+0] * m[3*4+1] + m[0*4+2] * m[1*4+3] * m[2*4+0] * m[3*4+1] +
    m[0*4+3] * m[1*4+0] * m[2*4+2] * m[3*4+1] - m[0*4+0] * m[1*4+3] * m[2*4+2] * m[3*4+1] -
    m[0*4+2] * m[1*4+0] * m[2*4+3] * m[3*4+1] + m[0*4+0] * m[1*4+2] * m[2*4+3] * m[3*4+1] +
    m[0*4+3] * m[1*4+1] * m[2*4+0] * m[3*4+2] - m[0*4+1] * m[1*4+3] * m[2*4+0] * m[3*4+2] -
    m[0*4+3] * m[1*4+0] * m[2*4+1] * m[3*4+2] + m[0*4+0] * m[1*4+3] * m[2*4+1] * m[3*4+2] +
    m[0*4+1] * m[1*4+0] * m[2*4+3] * m[3*4+2] - m[0*4+0] * m[1*4+1] * m[2*4+3] * m[3*4+2] -
    m[0*4+2] * m[1*4+1] * m[2*4+0] * m[3*4+3] + m[0*4+1] * m[1*4+2] * m[2*4+0] * m[3*4+3] +
    m[0*4+2] * m[1*4+0] * m[2*4+1] * m[3*4+3] - m[0*4+0] * m[1*4+2] * m[2*4+1] * m[3*4+3] -
    m[0*4+1] * m[1*4+0] * m[2*4+2] * m[3*4+3] + m[0*4+0] * m[1*4+1] * m[2*4+2] * m[3*4+3];
}

//m needs to be 3by3 row major, computes m.dot(p)
template<class pointT>
point<3> dot3(floatT* m, point<3> p) {
  point<3> r;
  r[0] = m[0]*p[0] + m[1]*p[1] + m[2]*p[2];
  r[1] = m[3]*p[0] + m[4]*p[1] + m[5]*p[2];
  r[2] = m[6]*p[0] + m[7]*p[1] + m[8]*p[2];
  return r;
}

//m needs to be 4by4 row major, computes m.dot(p)
template<class pointT>
point<4> dot4(floatT* m, point<4> p) {
  point<4> r;
  r[0] = m[0]*p[0] + m[1]*p[1] + m[2]*p[2] + m[3]*p[3];
  r[1] = m[4]*p[0] + m[5]*p[1] + m[6]*p[2] + m[7]*p[3];
  r[2] = m[8]*p[0] + m[9]*p[1] + m[10]*p[2] + m[11]*p[3];
  r[3] = m[12]*p[0] + m[13]*p[1] + m[14]*p[2] + m[15]*p[3];
  return r;
}

template<class pointT>
point<3> crossProduct(point<3> a, point<3> b) {
  point<3> r;
  r.updateCoordinate(0, a[1]*b[2]-a[2]*b[1]);
  r.updateCoordinate(1, a[2]*b[0]-a[0]*b[2]);
  r.updateCoordinate(2, a[0]*b[1]-a[1]*b[0]);
  return r;
}

// A class for sphere
class sphere {
public:
  typedef double floatT;
  typedef point<3> pointT;
  pointT ct;
  floatT rad;

  inline floatT radius() {return rad;}
  inline pointT center() {return ct;}
  floatT& operator[] (int i) {return ct[i];}

  sphere(pointT a, pointT b, pointT c, pointT d) {
    //reference: http://www.ambrsoft.com/TrigoCalc/Sphere/Spher3D_.htm
    auto x1 = a[0]; auto y1 = a[1]; auto z1 = a[2];
    auto x2 = b[0]; auto y2 = b[1]; auto z2 = b[2];
    auto x3 = c[0]; auto y3 = c[1]; auto z3 = c[2];
    auto x4 = d[0]; auto y4 = d[1]; auto z4 = d[2];
    floatT t1 = -(x1*x1+y1*y1+z1*z1);
    floatT t2 = -(x2*x2+y2*y2+z2*z2);
    floatT t3 = -(x3*x3+y3*y3+z3*z3);
    floatT t4 = -(x4*x4+y4*y4+z4*z4);
    floatT M1[16] = {x1,y1,z1,1,x2,y2,z2,1,x3,y3,z3,1,x4,y4,z4,1};
    floatT T = determinant4by4(M1);
    floatT M2[16] = {t1,y1,z1,1,t2,y2,z2,1,t3,y3,z3,1,t4,y4,z4,1};
    floatT D = determinant4by4(M2)/T;
    floatT M3[16] = {x1,t1,z1,1,x2,t2,z2,1,x3,t3,z3,1,x4,t4,z4,1};
    floatT E = determinant4by4(M3)/T;
    floatT M4[16] = {x1,y1,t1,1,x2,y2,t2,1,x3,y3,t3,1,x4,y4,t4,1};
    floatT F = determinant4by4(M4)/T;
    floatT M5[16] = {x1,y1,z1,t1,x2,y2,z2,t2,x3,y3,z3,t3,x4,y4,z4,t4};
    floatT G = determinant4by4(M5)/T;
    ct.x[0] = -D/2;
    ct.x[1] = -E/2;
    ct.x[2] = -F/2;
    rad = 0.5*sqrt(D*D+E*E+F*F-4*G);
  }

  sphere(pointT a, pointT b, pointT c) {
    //change basis
    auto v1 = (b-a).normalize();
    auto v3 = crossProduct<pointT>(c-a,b-a).normalize();
    auto v2 = crossProduct<pointT>(v1,v3);
    floatT T[9];
    T[0] = v1[0]; T[1] = v2[0]; T[2] = v3[0];
    T[3] = v1[1]; T[4] = v2[1]; T[5] = v3[1];
    T[6] = v1[2]; T[7] = v2[2]; T[8] = v3[2];
    floatT Ti[9]; for (int i=0; i<9; ++i) Ti[i] = T[i];
    inverse3by3(Ti);
    auto aa = dot3<pointT>(Ti, a);
    auto bb = dot3<pointT>(Ti, b);
    auto cc = dot3<pointT>(Ti, c);
    //2d circle
    auto x1 = aa[0]; auto y1 = aa[1];
    auto x2 = bb[0]; auto y2 = bb[1];
    auto x3 = cc[0]; auto y3 = cc[1];
    floatT A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2;
    floatT B = (x1*x1+y1*y1)*(y3-y2) + (x2*x2+y2*y2)*(y1-y3) + (x3*x3+y3*y3)*(y2-y1);
    floatT C = (x1*x1+y1*y1)*(x2-x3) + (x2*x2+y2*y2)*(x3-x1) + (x3*x3+y3*y3)*(x1-x2);
    floatT D = (x1*x1+y1*y1)*(x3*y2-x2*y3) + (x2*x2+y2*y2)*(x1*y3-x3*y1) + (x3*x3+y3*y3)*(x2*y1-x1*y2);
    ct.x[0] = -B/(2*A);
    ct.x[1] = -C/(2*A);
    rad = sqrt((B*B+C*C-4*A*D)/(4*A*A));
    //change back
    ct.x[2] = aa[2];
    ct = dot3<pointT>(T, ct);
  }

  sphere(pointT a, pointT b) {
    ct = a.average(b);
    rad = a.pointDist(b)/2;
  }

  sphere(pointT centerr, floatT radd): ct(centerr), rad(radd) {}

  //empty constructor
  sphere(): ct(pointT()), rad(-1) {}

  bool isEmpty() {return rad < 0;}

  inline bool contain(pointT p) {
    return p.pointDist(ct) <= rad*1.000001;}//todo, sqrt optimize
};

// A class for circle
class circle {
public:
  typedef double floatT;
  typedef point<2> pointT;
  pointT ct;
  floatT rad;

  inline floatT radius() {return rad;}
  inline pointT center() {return ct;}
  floatT& operator[] (int i) {return ct[i];}

  circle(pointT a, pointT b, pointT c) {
    //reference: http://www.ambrsoft.com/trigocalc/circle3d.htm
    auto x1 = a[0]; auto y1 = a[1];
    auto x2 = b[0]; auto y2 = b[1];
    auto x3 = c[0]; auto y3 = c[1];
    floatT A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2;
    floatT B = (x1*x1+y1*y1)*(y3-y2) + (x2*x2+y2*y2)*(y1-y3) + (x3*x3+y3*y3)*(y2-y1);
    floatT C = (x1*x1+y1*y1)*(x2-x3) + (x2*x2+y2*y2)*(x3-x1) + (x3*x3+y3*y3)*(x1-x2);
    floatT D = (x1*x1+y1*y1)*(x3*y2-x2*y3) + (x2*x2+y2*y2)*(x1*y3-x3*y1) + (x3*x3+y3*y3)*(x2*y1-x1*y2);
    ct.x[0] = -B/(2*A);
    ct.x[1] = -C/(2*A);
    rad = sqrt((B*B+C*C-4*A*D)/(4*A*A));
  }

  circle(pointT a, pointT b) {
    ct = a.average(b);
    rad = a.pointDist(b)/2;
  }

  circle(pointT centerr, floatT radd): ct(centerr), rad(radd) {}

  //empty constructor
  circle(): ct(pointT()), rad(-1) {}

  bool isEmpty() {return rad < 0;}

  inline bool contain(pointT p) {
    return p.pointDist(ct) <= rad*1.000001;}//todo, sqrt optimize
};

inline vect3d onParabola(vect2d v) {
  return vect3d(v.x(), v.y(), v.x()*v.x() + v.y()*v.y());}

// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c. 
// Projects a, b, c onto a parabola centered with d at the origin
//   and does a plane side test (tet volume > 0 test)
inline bool inCircle(point2d a, point2d b, 
                     point2d c, point2d d) {
  vect3d ad = onParabola(a-d);
  vect3d bd = onParabola(b-d);
  vect3d cd = onParabola(c-d);
  return (ad.cross(bd)).dot(cd) > 0.0;
}

// returns a number between -1 and 1, such that -1 is out at infinity,
// positive numbers are on the inside, and 0 is at the boundary
inline double inCircleNormalized(point2d a, point2d b,
				 point2d c, point2d d) {
  vect3d ad = onParabola(a-d);
  vect3d bd = onParabola(b-d);
  vect3d cd = onParabola(c-d);
  return (ad.cross(bd)).dot(cd)/(ad.Length()*bd.Length()*cd.Length());
}

// *************************************************************
//    TRIANGLES
// *************************************************************

struct triangle {
  int C[3];
  triangle(int p1, int p2, int p3) {
    C[0] = p1; C[1] = p2; C[2] = p3;
  }
};

template <class point>
struct triangles {
  intT numPoints;
  intT numTriangles;
  point* P;
  triangle* T;
  triangles() {}
  void del() {free(P); free(T);}
  triangles(intT np, intT nt, point* _P, triangle* _T) 
    : numPoints(np), numTriangles(nt), P(_P), T(_T) {}
};

template <class pointT>
struct ray {
  typedef typename pointT::vectT vectT;
  pointT o;
  vectT d;
  ray(pointT _o, vectT _d) : o(_o), d(_d) {}
  ray() {}
};

inline double angle(point2d a, point2d b, point2d c) {
  vect2d ba = (b-a);
  vect2d ca = (c-a);
  double lba = ba.Length();
  double lca = ca.Length();
  double pi = 3.14159;
  return 180/pi*acos(ba.dot(ca)/(lba*lca));
}

inline double minAngle(point2d a, point2d b, point2d c) {
  double aAngle = angle(a,b,c);
  double cAngle = angle(c,b,a);
  return min(180-aAngle-cAngle, min(aAngle,cAngle));
}

inline double minAngleCheck(point2d a, point2d b, point2d c, double angle) {
  vect2d ba = (b-a);
  vect2d ca = (c-a);
  vect2d cb = (c-b);
  double lba = ba.Length();
  double lca = ca.Length();
  double lcb = cb.Length();
  double pi = 3.14159;
  double co = cos(angle*pi/180.);
  return (ba.dot(ca)/(lba*lca) > co || ca.dot(cb)/(lca*lcb) > co || 
	  -ba.dot(cb)/(lba*lcb) > co);
}

inline point2d triangleCircumcenter(point2d a, point2d b, point2d c) {
  vect2d v1 = b-a;
  vect2d v2 = c-a;
  vect2d v11 = v1 * v2.dot(v2);
  vect2d v22 = v2 * v1.dot(v1);
  return  a + vect2d(v22.y() - v11.y(), v11.x() - v22.x())/(2.0 * v1.cross(v2));
}

#endif // _BENCH_GEOM_INCLUDED
