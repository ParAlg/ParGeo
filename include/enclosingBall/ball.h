#pragma once

#include "pargeo/point.h"
#include "parlay/sequence.h"

#define ROBUST_BALL // shifting primitive magnitude to reduce numerical error

namespace pargeo {
  namespace seb {
    template <int dim> class ball;
  }
}

template <int dim>
class pargeo::seb::ball {
  using pointT = pargeo::point<dim>;
  using floatT = typename pointT::floatT;

  int d;//number of supporting points
  floatT AInv[dim*dim];
  pointT P[dim+1];
  floatT Q[dim];
  floatT La[dim];
#ifdef ROBUST_BALL
  pointT offset;
#endif

  pointT c;
  floatT r;

  //takes row-major 2by2 matrix m
  inline floatT determinant2by2(floatT* m) {
    return m[0]*m[3] - m[1]*m[2];
  }

  //takes row-major 2by2 matrix m, output replaces m
  inline void inverse2by2(floatT* m) {
    auto coeff = 1/determinant2by2(m);
    floatT mInv[4];
    mInv[0] = m[3]*coeff;
    mInv[1] = -m[1]*coeff;
    mInv[2] = -m[2]*coeff;
    mInv[3] = m[0]*coeff;
    m[0] = mInv[0]; m[1] = mInv[1]; m[2] = mInv[2]; m[3] = mInv[3];
  }

  void recompute() {
    for(int i=0; i<dim; ++i) c[i] = 0;
    for(int i=1; i<d; ++i) {
      c = c + (P[i] - P[0]).mult(La[i-1]);
    }
    r = sqrt(c.dot(c));
    c = c+P[0];
  }

  void twoPointConstruct() {
    //c = P[0].average(P[1]);
    c = (P[0] + P[1])/2;
    r = P[0].dist(P[1])/2;
  }

  void threePointConstruct() {
    floatT* A = AInv;
    for(int i=1; i<d; ++i) {
      auto Q1 = P[i] - P[0];
      Q[i-1] = Q1.dot(Q1);
      for (int j=1; j<d; ++j) {
        auto Q2 = P[j] - P[0];
        A[(i-1)*(d-1)+(j-1)] = 2*Q1.dot(Q2);
      }
    }

    //we knew d=3, and A is 2x2
    inverse2by2(A);//A=AInv

    //Lambda=AInv.dot(Q)
    for(int i=0; i<(d-1); ++i) {
      La[i] = 0;
      for(int j=0; j<(d-1); ++j) {
        La[i] += AInv[i*(d-1)+j]*Q[j];}
    }
    recompute();
  }

public:
#ifdef ROBUST_BALL
  inline pointT center() {return c+offset;}
  inline bool contain(pointT p) {
    auto pp = p-offset;
    //return pp.dist(c) <= radius()+1e-9;}
    return pp.dist(c) <= radius()+p.eps;}
#else
  inline pointT center() {return c;}
  inline bool contain(pointT p) {
    //return p.dist(c) <= radius()+1e-9;}
    return p.dist(c) <= radius()+p.eps;}
#endif
  inline pointT* support() {return P;}
  inline floatT radius() {return r;}
  inline int size() {return d;}
  inline bool isEmpty() {return size() <= 0;}

  ball(): d(0) {}

  //ball(pointT* PP, int dd): d(dd) {
  ball(parlay::slice<pointT*,pointT*> PP): d(PP.size()) {
#ifdef ROBUST_BALL
    for(int i=0; i<dim; ++i) offset[i] = 0;
    for(int i=0; i<d; ++i) {
      P[i] = PP[i];
      offset = offset + PP[i];
    }
    offset = offset/d;
    for(int i=0; i<d; ++i) P[i] = P[i]-offset;//offset supports to a safer zone
#else
    for(int i=0; i<d; ++i) P[i] = PP[i];
#endif
    if (d <= 1) {
      throw std::runtime_error("cannot construct ball on <=1 point");
    } else if (d == 2) twoPointConstruct();
    else if (d == 3) threePointConstruct();
    else if (d > 3) {
      //std::cout << "d>3 construction\n";
      d = 3;
      threePointConstruct();
      for (int i=3; i<PP.size(); ++i) grow(P[i]);
    } else {
      throw std::runtime_error("ball wrong dimension, abort");
    }
#ifdef ROBUST_BALL
    for(int i=0; i<d; ++i) P[i] = P[i]+offset;//restore supports
#endif
  }

  void grow(pointT q) {
    //std::cout << "grow " << q << "\n";
    if (d+1 > dim+1) {
      std::cout << "d+1 = " << d+1 << "\n";
      throw std::runtime_error("ball max points exceeded");
    }

    P[d] = q;
    d += 1;
    if (d == 3) {
      return threePointConstruct();
    }

    //expand inverse of A, recompute center and radius
    floatT mu[d-2];
    floatT tmp[d-2];//the Q vector
    point<dim> Qm = P[d-1] - P[0];

    for(int j=1; j<d-1; ++j) {
      tmp[j-1] = 2*(P[j]-P[0]).dot(Qm);}

    //mu=AInv.dot([2*Q1.dot(Qm) ... 2*Qm-1.dot(Qm)])
    for(int i=0; i<(d-2); ++i) {
      mu[i] = 0;
      for(int j=0; j<(d-2); ++j) {
        mu[i] += AInv[i*(d-2)+j]*tmp[j];}
    }

    //z=
    floatT tmp2 = 0;
    for(int i=0; i<(d-2); ++i) tmp2 += tmp[i]*mu[i];
    floatT z = 2*Qm.dot(Qm) - tmp2;

    floatT newAInv[dim*dim];
    for (int i=0; i<d-2; ++i) {
      for (int j=0; j<d-2; ++j) {
        newAInv[i*(d-1)+j] = AInv[i*(d-2)+j] + mu[i]*mu[j]/z;}
    }
    for (int i=0; i<d-2; ++i) {
      newAInv[i*(d-1)+d-2] = -mu[i]/z;
      newAInv[(d-2)*(d-1)+i] = -mu[i]/z;
    }
    newAInv[(d-2)*(d-1)+(d-2)] = 1/z;
    for(int i=0; i<(d-1)*(d-1); ++i) {
      AInv[i] = newAInv[i];}

    //Lambda=AInv.dot(Q)
    Q[d-2] = (P[d-1]-P[0]).dot(P[d-1]-P[0]);//update Q
    for(int i=0; i<(d-1); ++i) {
      La[i] = 0;
      for(int j=0; j<(d-1); ++j) {
        La[i] += AInv[i*(d-1)+j]*Q[j];}
    }
    recompute();
  }

};
