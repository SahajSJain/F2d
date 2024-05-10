#pragma once
// #include <eigen3/Eigen/Dense>
#include "../eigen3/Eigen/Dense"

namespace CLASSES {
class SolverCoeffs { // SolverCoeffs for Ax=b//nodewise
public:
  Eigen::ArrayXXd E, W, S, N, P, Pinv; //
  SolverCoeffs() {                     // Empty Constructor does nothing
  }
  void assign(const int &NX, const int &NY) {
    E = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    W = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    N = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    S = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    P = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    Pinv = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
  }
  void
  operator=(const SolverCoeffs &B) { // overload '=' operator. Allows to do A=B;
    E = B.E;
    W = B.W;
    N = B.N;
    S = B.S;
    P = B.P;
    Pinv = B.Pinv;
  }
};
class PrimFlowField { // Primitive variables for flow fiels
public:
  Eigen::ArrayXXd u, us, v, vs, p, p_pert, ci, cj, u_i, u_j, cis, cjs; //
  Eigen::ArrayXXd om, omx,
      omy; // I GUESS Vorticity is not primitive // The idea was originally to
           // have a secondary flow field with values like phase and vorticityt
  PrimFlowField() { // Empty Constructor does nothing
  }
  void operator=(
      const PrimFlowField &P) { // overload '=' operator. Allows to do A=B;
    u = P.u;
    us = P.us;
    v = P.v;
    vs = P.vs;
    p = P.p;
    p_pert = P.p_pert;
    ci = P.ci;
    cj = P.cj;
    u_i = P.u_i;
    u_j = P.u_j;
    cis = P.cis;
    cjs = P.cjs;
    om = P.om;
    omx = P.omx;
    omy = P.omy;
  }
  void InitializeCenters(int &NX, int &NY, double &uinit, double &vinit,
                         double &pert) {
    u = Eigen::ArrayXXd::Constant(NX, NY, uinit) +
        pert * Eigen::ArrayXXd::Random(NX, NY);
    v = Eigen::ArrayXXd::Constant(NX, NY, vinit) +
        pert * Eigen::ArrayXXd::Random(NX, NY);
    us = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    vs = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    p = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    p_pert = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    omx = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    omy = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
    om = Eigen::ArrayXXd::Constant(NX, NY, 0.0);
  }
  void InitializeFaces(int &NX, int &NY, double &uinit, double &vinit) {
    u_i = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    u_j = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    ci = Eigen::ArrayXXd::Constant(NX, NY, uinit);
    cj = Eigen::ArrayXXd::Constant(NX, NY, vinit);
    cis = Eigen::ArrayXXd::Constant(NX, NY, 0);
    cjs = Eigen::ArrayXXd::Constant(NX, NY, 0);
  }
};

class FractionalStepArrays { //
public:
  Eigen::ArrayXXd gradpx, gradpy;
  Eigen::ArrayXXd uD1, uC1, uC0;
  Eigen::ArrayXXd vD1, vC1, vC0;
  Eigen::ArrayXXd RHSu, RHSv;
  Eigen::ArrayXXd div_s, div_new, f;
  Eigen::ArrayXXd gradpx_f, gradpy_f; //
  Eigen::MatrixXd Res1, Res2;
  FractionalStepArrays() { // Empty Constructor does nothing
  }
  void InitializeCenters(int &NX, int &NY) {
    gradpx = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    gradpy = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    uD1 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    uC1 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    uC0 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    vD1 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    vC1 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    vC0 = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    RHSu = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    RHSv = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    Res1 = Eigen::MatrixXd::Constant(NX, NY, 0.);
    Res2 = Eigen::MatrixXd::Constant(NX, NY, 0.);

    div_s = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    div_new = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    f = Eigen::ArrayXXd::Constant(NX, NY, 0.);
  }
  void InitializeFaces(int &NX, int &NY) {
    gradpx_f = Eigen::ArrayXXd::Constant(NX, NY, 0.);
    gradpy_f = Eigen::ArrayXXd::Constant(NX, NY, 0.);
  }
};

class NeighbourField { // SolverCoeffs for Ax=b//nodewise
public:
  Eigen::ArrayXXd E, W, S, N, P; //
  NeighbourField() {             // Empty Constructor does nothing
  }
  void assign(const int &NX, const int &NY) {
    P = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    E = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    W = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    N = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
    S = Eigen::ArrayXXd::Constant(NX, NY, 1.0);
  }
  void operator=(
      const NeighbourField &B) { // overload '=' operator. Allows to do A=B;
    P = B.P;
    E = B.E;
    W = B.W;
    N = B.N;
    S = B.S;
  }
};

class BlankArray { // SolverCoeffs for Ax=b//nodewise
public:
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ibnode, blank, E, W, S,
      N; //

  BlankArray() { // Empty Constructor does nothing
  }
  void assign(const int &NX, const int &NY) {
    blank =
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 0);
    E = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 1);
    W = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 1);
    N = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 1);
    S = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 1);
    ibnode =
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(NX, NY, 0);
  }
  void
  operator=(const BlankArray &B) { // overload '=' operator. Allows to do A=B;
    E = B.E;
    W = B.W;
    N = B.N;
    S = B.S;
    blank = B.blank;
    ibnode = B.blank;
  }
};

class BoundaryArray { // SolverCoeffs for Ax=b//nodewise
public:
  Eigen::VectorXi i, j; // Inverse of dels
  int Nb;
  Eigen::VectorXd xc, yc, fx, fy; // Inverse of dels
  double forcex, forcey, moment;
  double xr, yr;
  int ih, jh;
  BoundaryArray() { // Empty Constructor does nothing
  }
  void assign(const int &Nx) {
    Nb = Nx;
    i = Eigen::VectorXi::Zero(Nb);
    j = Eigen::VectorXi::Zero(Nb);
    xc = Eigen::VectorXd::Zero(Nb);
    yc = Eigen::VectorXd::Zero(Nb);
    fx = Eigen::VectorXd::Zero(Nb);
    fy = Eigen::VectorXd::Zero(Nb);
  }
  void operator=(
      const BoundaryArray &B) { // overload '=' operator. Allows to do A=B;
    Nb = B.Nb;
    i = B.i;
    j = B.j;
  }
  void totalforce(double &xp, double &yp) {
    forcex = 0.0;
    forcey = 0.0;
    moment = 0.0;
    int issmh = 0;
#pragma omp parallel for private(issmh) reduction(+ : forcex, forcey, moment)  \
    shared(fx, fy, xc, yc , xp , yp)
    for (issmh = 0; issmh < Nb; issmh++) {
      xr=xc(issmh)-xp; //relative x and y position
      yr=yc(issmh)-yp;
      forcex = forcex + fx(issmh);
      forcey = forcey + fy(issmh);
      moment = moment + fy(issmh)*xr-fx(issmh)*yr;
    }
  }
};
} // namespace CLASSES