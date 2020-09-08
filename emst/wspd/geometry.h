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
  vectT operator-(pointT op2) {
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

// A class for sphere
template <int _dim> class sphere {
public:
  typedef double floatT;
  typedef point<_dim> pointT;
  pointT center;
  floatT radius;

  //_dim=2 circle constructor
  sphere(pointT a, pointT b, pointT c) {
    if (_dim != 2) {
      cout << "error, this constructor for sphere only works for dim=2, abort" << endl;}
    auto x1 = a[0]; auto y1 = a[1];
    auto x2 = b[0]; auto y2 = b[1];
    auto x3 = c[0]; auto y3 = c[1];
    floatT A = x1*(y2-y3) - y1*(x2-x3) + x2*y3 - x3*y2;
    floatT B = (x1*x1+y1*y1)*(y3-y2) + (x2*x2+y2*y2)*(y1-y3) + (x3*x3+y3*y3)*(y2-y1);
    floatT C = (x1*x1+y1*y1)*(x2-x3) + (x2*x2+y2*y2)*(x3-x1) + (x3*x3+y3*y3)*(x1-x2);
    floatT D = (x1*x1+y1*y1)*(x3*y2-x2*y3) + (x2*x2+y2*y2)*(x1*y3-x3*y1) + (x3*x3+y3*y3)*(x2*y1-x1*y2);
    center.x[0] = -B/(2*A);
    center.x[1] = -C/(2*A);
    radius = sqrt((B*B+C*C-4*A*D)/(4*A*A));
  }

  //generic constructor 1
  sphere(pointT centerr, floatT radiuss): center(centerr), radius(radiuss) {}

  //generic constructor 2
  sphere(pointT a, pointT b) {
    center = a.average(b);
    radius = a.pointDist(b)/2;
  }

  //empty constructor
  sphere(): center(pointT()), radius(-1) {}

  inline bool contain(pointT p) {
    return p.pointDist(center) <= radius*1.000001;}//todo, sqrt optimize
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

