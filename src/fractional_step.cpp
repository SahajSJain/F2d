#include <eigen3/Eigen/Dense>
#include <iostream>
// #include <string> //Wheres Wallace????make
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"

// #include "./headers/FIELDVARS.h"
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FLOWVARS.h"
#include "headers/INPUT.h"


// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
// using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;

void fractional_step() {

  // start iterations, non VK
  Older = Old;
  Old = New;
  New.us = New.u; /// initializa
  New.vs = New.v; /// initialize
  // Calculate pressure gradient
  //    GradPFace(Old.p, FSA.gradpx_f, FSA.gradpy_f);
  GradPCenter(Old.p, FSA.gradpx, FSA.gradpy);

  // Assemble RHS for us
  DiffCalc(FSA.uD1, Old.u, DEq); // Calculate D(u_n)
  ConvCalc(FSA.uC1, Old.u, Old.u_i, Old.u_j, Old.ci,
           Old.cj); // Calculates NL(u_n)
  ConvCalc(FSA.uC0, Older.u, Older.u_i, Older.u_j, Older.ci,
           Older.cj); // Calculates NL(u_n-1)

  FSA.RHSu = -1.5 * FSA.uC1 + 0.5 * FSA.uC0 + 0.5 * Reinv * FSA.uD1 +
             Old.u * dtinv - vkflag * FSA.gradpx;
  FSA.RHSu += force_x; // apply force in x firection
  // Assemble RHS for vs
  DiffCalc(FSA.vD1, Old.v, DEq); // Calculate D(u_n)
  ConvCalc(FSA.vC1, Old.v, Old.u_i, Old.u_j, Old.ci,
           Old.cj); // Calculates NL(u_n)
  ConvCalc(FSA.vC0, Older.v, Older.u_i, Older.u_j, Older.ci,
           Older.cj); // Calculates NL(u_n-1)

  FSA.RHSv = -1.5 * FSA.vC1 + 0.5 * FSA.vC0 + 0.5 * Reinv * FSA.vD1 +
             Old.v * dtinv - vkflag * FSA.gradpy;
  FSA.RHSu += force_y; // apply force in x firection

  // cout<<Old.v;

  // Solve for us and vs

  int kmom = 0;
  double errormom = 0;
  VelocityPredictor(New.us, New.vs, FSA.RHSu, FSA.RHSv, ADC, kmom, errormom);
  //   cout<<Old.u;
  if (printflag == 0)
    outputout << "\n"
              << "\n"
              << "Mom Iteration count:" << kmom << "\n";
  if (printflag == 0)
    outputout << "Mom Error:" << errormom << "\n"
              << "\n";
  // cout << "\n"<<New.ci<<"\n"<<"\n"<<"\n";
  HalfCalc(New.us, New.cis, New.u_j); // Calculate cis
  HalfCalc(New.vs, New.u_i, New.cjs); // Calculate cjs
  // apply bc on intermediate face fluxes
  FaceFluxBC(New.cis, New.cjs);

  // calculate divergence of (cis,cjs) field
  MassCalc(FSA.div_s, New.cis,
           New.cjs); // div_s is the intermeediate divergence
  double massold = FSA.div_s.matrix().norm() / (Nx * Ny);
  if (printflag == 0)
    outputout << "\n"
              << "\n"
              << "Intermediate Mass total: " << massold << "\n";

  FSA.f = FSA.div_s * dtinv;
  // Solve for pressure
  int kPsolver = 0;
  double Perror = 0;
  PointSOR(New.p_pert, FSA.f, DEq, kPsolver, Perror);
  // Calculate pressure gradient
  GradPFace(New.p_pert, FSA.gradpx_f, FSA.gradpy_f);
  GradPCenter(New.p_pert, FSA.gradpx, FSA.gradpy);
  // Correct velocity and flux field
  New.u = New.us - dt * FSA.gradpx;
  New.v = New.vs - dt * FSA.gradpy;
  New.ci = New.cis - dt * FSA.gradpx_f;
  New.cj = New.cjs - dt * FSA.gradpy_f;
  // apply boundary conditions for corrected velocity field
  // VelBC(New.u, New.v);
  // FaceFluxBC(New.ci, New.cj);
  MassCalc(FSA.div_new, New.ci, New.cj); // div_new is the corrected divergence
  // cout<<"\n"<<"\n"<<"Intermediate divergence
  // field"<<"\n"<<"\n"<<FSA.div_s<<"\n"; cout<<"\n"<<"\n"<<"Final divergence
  // field"<<"\n"<<"\n"<<FSA.div_new<<"\n";
  double massnew = FSA.div_new.matrix().norm() / (Nx * Ny);
  // Update pressure
  New.p = New.p_pert + vkflag * Old.p;
  if (printflag == 0)
    outputout << "\n"
              << "\n"
              << "Pressure Iteration count:" << kPsolver << "\n";
  if (printflag == 0)
    outputout << "Pressure Error:" << Perror << "\n";
  if (printflag == 0)
    outputout << "Mass total after correction = " << massnew << "\n";

  //    if(printflag==0) outputout << "Mom Error:" << errormom << "\n"<<"\n";

  //    cout << "Old Mass total: " << massold << "\n";
  //    cout << "New Mass total: " << massnew << "\n";
  //    cout << "Mom Iteration count:" << kmom << "\n";
  //    cout << "Pressure Iteration count:" << kPsolver << "\n";
  //
  //    cout << "Mom Error:" << errormom << "\n";
  //    cout << "Pressure Error:" << Perror << "\n";
  double velnorm =
      (New.u - Old.u).matrix().norm() + (New.v - Old.v).matrix().norm();
  velnorm = velnorm / (Nx * Ny);
  PresBC(New.p);
  if (debugflag == 1)
    cout << "\n"
         << timestep << "\t" << soltime << "\t" << kPsolver << "\t" << massold
         << "\t" << massnew << "\t" << velnorm << "\t" << fluxintegral << "\t"
         << flowintegralpost << "\t" << outflowarea << "\n";
//   if (debugflag == 1)
//     cout << "\n"
//          << inflowintegral << "\t" << outflowintegral << "\t" 
//          << "\t" << inflowintegralpost << "\t" << outflowintegralpost << "\n";
  // tecplot_printer();
}